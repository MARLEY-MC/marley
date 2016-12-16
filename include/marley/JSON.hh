// Based on https://github.com/lefticus/SimpleJSON
#pragma once

// standard library includes
#include <cctype>
#include <cmath>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <iostream>
#include <istream>
#include <ostream>
#include <limits>
#include <sstream>
#include <map>
#include <string>
#include <type_traits>

// MARLEY includes
#include "marley/Logger.hh"

namespace marley {

  namespace {
    std::string json_escape( const std::string& str ) {
      std::string output;
      for( unsigned i = 0; i < str.length(); ++i )
      switch( str[i] ) {
        case '\"': output += "\\\""; break;
        case '\\': output += "\\\\"; break;
        case '\b': output += "\\b";  break;
        case '\f': output += "\\f";  break;
        case '\n': output += "\\n";  break;
        case '\r': output += "\\r";  break;
        case '\t': output += "\\t";  break;
        default  : output += str[i]; break;
      }
      return output;
    }
  }

  class JSON
  {
    union Data {
      Data(double d) : float_(d) {}
      Data(long l) : integer_(l) {}
      Data(bool b) : boolean_(b) {}
      Data(const std::string& s) : string_(new std::string(s)) {}
      Data() : integer_(0) {}

      std::deque<JSON>* list_;
      std::map<std::string, JSON>* map_;
      std::string* string_;
      double float_;
      long integer_;
      bool boolean_;
    } data_;

    public:
      enum class DataType {
        Null,
        Object,
        Array,
        String,
        Floating,
        Integral,
        Boolean
      };

      template <typename Container> class JSONWrapper {

        private:
          Container* object;

        public:
          JSONWrapper(Container* val) : object(val) {}
          JSONWrapper(std::nullptr_t) : object(nullptr) {}

          typename Container::iterator begin() {
            return object ? object->begin() : typename Container::iterator();
          }
          typename Container::iterator end() {
            return object ? object->end() : typename Container::iterator();
          }
          typename Container::const_iterator begin() const {
            return object ? object->begin() : typename Container::iterator();
          }
          typename Container::const_iterator end() const {
            return object ? object->end() : typename Container::iterator();
          }
      };

      template <typename Container>
      class JSONConstWrapper {

        private:
          const Container* object;

        public:
          JSONConstWrapper(const Container* val) : object(val) {}
          JSONConstWrapper(std::nullptr_t) : object(nullptr) {}

          typename Container::const_iterator begin() const
            { return object ? object->begin()
              : typename Container::const_iterator(); }
          typename Container::const_iterator end() const
            { return object ? object->end()
              : typename Container::const_iterator(); }
      };

      JSON() : data_(), type_(DataType::Null) {}

      JSON(std::initializer_list<JSON> list) : JSON()
      {
        set_type(DataType::Object);
        for(auto i = list.begin(), e = list.end(); i != e; ++i, ++i)
          operator[](i->to_string()) = *std::next(i);
      }

      JSON(JSON&& other) : data_(other.data_), type_(other.type_)
      { other.type_ = DataType::Null; other.data_.map_ = nullptr; }

      JSON& operator=(JSON&& other) {
        data_ = other.data_;
        type_ = other.type_;
        other.data_.map_ = nullptr;
        other.type_ = DataType::Null;
        return *this;
      }

      JSON(const JSON& other) {
        switch(other.type_) {
        case DataType::Object:
          data_.map_ = new std::map<std::string,JSON>(
            other.data_.map_->begin(), other.data_.map_->end());
          break;
        case DataType::Array:
          data_.list_ = new std::deque<JSON>(other.data_.list_->begin(),
            other.data_.list_->end());
          break;
        case DataType::String:
          data_.string_ = new std::string(*other.data_.string_);
          break;
        default:
          data_ = other.data_;
        }
        type_ = other.type_;
      }

      JSON& operator=(const JSON& other) {
        switch(other.type_) {
          case DataType::Object:
            data_.map_ = new std::map<std::string,JSON>(
              other.data_.map_->begin(), other.data_.map_->end());
            break;
          case DataType::Array:
            data_.list_ = new std::deque<JSON>( other.data_.list_->begin(),
              other.data_.list_->end());
            break;
          case DataType::String:
            data_.string_ = new std::string(*other.data_.string_);
            break;
          default:
            data_ = other.data_;
        }
        type_ = other.type_;
        return *this;
      }

      ~JSON() {
        switch(type_) {
          case DataType::Array:
            delete data_.list_;
            break;
          case DataType::Object:
            delete data_.map_;
            break;
          case DataType::String:
            delete data_.string_;
            break;
          default:;
        }
      }

      template <typename T> JSON(T b,
        typename std::enable_if<std::is_same<T,bool>::value>::type* = 0)
        : data_(b), type_(DataType::Boolean) {}

      template <typename T> JSON(T i,
        typename std::enable_if<std::is_integral<T>::value
        && !std::is_same<T,bool>::value>::type* = 0)
        : data_(static_cast<long>(i)), type_(DataType::Integral) {}

      template <typename T> JSON(T f,
        typename std::enable_if<std::is_floating_point<T>::value>::type* = 0)
        : data_(static_cast<double>(f)), type_(DataType::Floating) {}

      template <typename T> JSON(T s,
        typename std::enable_if<std::is_convertible<T,
        std::string>::value>::type* = 0) : data_(std::string(s)),
        type_(DataType::String) {}

      JSON(std::nullptr_t) : data_(), type_(DataType::Null) {}

      static JSON make(DataType type) {
        JSON ret;
        ret.set_type(type);
        return ret;
      }

      static JSON load(const std::string& s);
      static JSON load(std::istream& is);
      static JSON load_file(const std::string& s);

      template <typename T> void append(T arg) {
        set_type(DataType::Array);
        data_.list_->emplace_back(arg);
      }

      template <typename T, typename... U> void append(T arg, U... args) {
        append(arg); append(args...);
      }

      template <typename T>
        typename std::enable_if<std::is_same<T,bool>::value, JSON&>::type
        operator=(T b)
      {
          set_type(DataType::Boolean);
          data_.boolean_ = b;
          return *this;
      }

      template <typename T> typename std::enable_if<std::is_integral<T>::value
        && !std::is_same<T,bool>::value, JSON&>::type operator=(T i)
      {
        set_type( DataType::Integral );
        data_.integer_ = i;
        return *this;
      }

      template <typename T>
        typename std::enable_if<std::is_floating_point<T>::value, JSON&>::type
        operator=(T f)
      {
        set_type(DataType::Floating);
        data_.float_ = f;
        return *this;
      }

      template <typename T> typename std::enable_if<std::is_convertible<T,
        std::string>::value, JSON&>::type operator=(T s)
      {
        set_type(DataType::String);
        *data_.string_ = std::string(s);
        return *this;
      }

      JSON& operator[](const std::string& key) {
        set_type(DataType::Object);
        return data_.map_->operator[](key);
      }

      JSON& operator[](unsigned index) {
        set_type(DataType::Array);
        if (index >= data_.list_->size()) data_.list_->resize(index + 1);
        return data_.list_->operator[](index);
      }

      JSON& at(const std::string& key) {
        return operator[](key);
      }

      const JSON& at(const std::string &key) const {
        return data_.map_->at(key);
      }

      JSON& at(unsigned index) {
        return operator[](index);
      }

      const JSON& at(unsigned index) const {
        return data_.list_->at(index);
      }

      int length() const {
        if (type_ == DataType::Array) return data_.list_->size();
        else return -1;
      }

      bool has_key(const std::string& key) const {
        if (type_ == DataType::Object)
          return data_.map_->find( key ) != data_.map_->end();
        else return false;
      }

      int size() const {
        if (type_ == DataType::Object)
          return data_.map_->size();
        else if (type_ == DataType::Array)
          return data_.list_->size();
        else
          return -1;
      }

      inline DataType type() const { return type_; }

      /// Functions for getting primitives from the JSON object.
      inline bool is_null() const { return type_ == DataType::Null; }
      inline bool is_object() const { return type_ == DataType::Object; }
      inline bool is_array() const { return type_ == DataType::Array; }
      inline bool is_string() const { return type_ == DataType::String; }
      inline bool is_float() const { return type_ == DataType::Floating; }
      inline bool is_integer() const { return type_ == DataType::Integral; }
      inline bool is_bool() const { return type_ == DataType::Boolean; }

      std::string to_string() const {
        bool b;
        return to_string(b);
      }

      std::string to_string(bool& ok) const {
        ok = (type_ == DataType::String);
        return ok ? json_escape(*data_.string_) : std::string("");
      }

      std::string to_string_or_throw() const {
        bool ok;
        std::string result = to_string(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value to string");
        return result;
      }

      double to_double() const {
        bool b;
        return to_double(b);
      }

      double to_double(bool& ok) const {
        ok = (type_ == DataType::Floating);
        if (ok) return data_.float_;
        ok = (type_ == DataType::Integral);
        if (ok) return data_.integer_;
        return 0.;
      }

      double to_double_or_throw() const {
        bool ok;
        double result = to_double(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to double");
        return result;
      }

      long to_long() const {
        bool b;
        return to_long( b );
      }

      long to_long(bool& ok) const {
        ok = (type_ == DataType::Integral);
        return ok ? data_.integer_ : 0;
      }

      long to_long_or_throw() const {
        bool ok;
        double result = to_long(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to long");
        return result;
      }

      bool to_bool() const {
        bool b;
        return to_bool( b );
      }

      bool to_bool(bool& ok) const {
        ok = (type_ == DataType::Boolean);
        return ok ? data_.boolean_ : false;
      }

      bool to_bool_or_throw() const {
        bool ok;
        double result = to_bool(ok);
        if (!ok) throw marley::Error("Failed to convert JSON value '"
          + to_string() + "' to bool");
        return result;
      }

      JSONWrapper<std::map<std::string,JSON> > object_range() {
        if (type_ == DataType::Object)
          return JSONWrapper<std::map<std::string,JSON>>(data_.map_);
        else return JSONWrapper<std::map<std::string,JSON>>(nullptr);
      }

      JSONWrapper<std::deque<JSON> > array_range() {
        if (type_ == DataType::Array)
          return JSONWrapper<std::deque<JSON>>(data_.list_);
        else return JSONWrapper<std::deque<JSON>>(nullptr);
      }

      JSONConstWrapper<std::map<std::string,JSON> > object_range() const {
        if (type_ == DataType::Object)
          return JSONConstWrapper<std::map<std::string,JSON>>(data_.map_);
        else return JSONConstWrapper<std::map<std::string,JSON>>(nullptr);
      }


      JSONConstWrapper<std::deque<JSON>> array_range() const {
        if ( type_ == DataType::Array )
          return JSONConstWrapper<std::deque<JSON>>(data_.list_);
        else return JSONConstWrapper<std::deque<JSON>>(nullptr);
      }

      std::string dump_string(int depth = 1, std::string tab = "  ") const {
        std::string pad("");
        for(int i = 0; i < depth; ++i, pad += tab);

        switch( type_ ) {
          case DataType::Null:
            return "null";
          case DataType::Object: {
             std::string s = "{\n";
            bool skip = true;
            for( auto &p : *data_.map_ ) {
              if ( !skip ) s += ",\n";
              s += ( pad + "\"" + p.first + "\" : "
                + p.second.dump_string( depth + 1, tab ) );
              skip = false;
            }
            s += ( "\n" + pad.erase( 0, 2 ) + "}" ) ;
            return s;
          }
          case DataType::Array: {
             std::string s = "[";
            bool skip = true;
            for( auto &p : *data_.list_ ) {
              if ( !skip ) s += ", ";
              s += p.dump_string( depth + 1, tab );
              skip = false;
            }
            s += "]";
            return s;
          }
          case DataType::String:
            return "\"" + json_escape( *data_.string_ ) + "\"";
          case DataType::Floating:
            return std::to_string( data_.float_ );
          case DataType::Integral:
            return std::to_string( data_.integer_ );
          case DataType::Boolean:
            return data_.boolean_ ? "true" : "false";
          default:
            return "";
        }
      }

    private:

      void check_if_object(const std::string& key) {
        if (type_ != DataType::Object) throw marley::Error("Attempted"
          " to retrieve a value for the key '" + key + " from a JSON primitive"
          " that is not an object");
      }

      void set_type(DataType type) {
        if (type == type_) return;

        switch(type_) {
          case DataType::Object:
            delete data_.map_;
            break;
          case DataType::Array:
            delete data_.list_;
            break;
          case DataType::String:
            delete data_.string_;
            break;
          default:;
        }

        switch(type) {
          case DataType::Null:
            data_.map_ = nullptr;
            break;
          case DataType::Object:
            data_.map_ = new std::map<std::string, JSON>();
            break;
          case DataType::Array:
            data_.list_ = new std::deque<JSON>();
            break;
          case DataType::String:
            data_.string_ = new std::string();
            break;
          case DataType::Floating:
            data_.float_ = 0.;
            break;
          case DataType::Integral:
            data_.integer_ = 0;
            break;
          case DataType::Boolean:
            data_.boolean_ = false;
            break;
        }

        type_ = type;
      }

    public:

      // Attempts to get a floating point number from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      double get_double(const std::string& key) {
        check_if_object(key);
        if (has_key(key)) return this->at(key).to_double_or_throw();
        else throw marley::Error("Missing JSON key '" + key + '\'');
        return 0.;
      }

      // Attempts to get a floating point number from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      double get_double(const std::string& key, double default_value) {
        check_if_object(key);
        if (!has_key(key)) return default_value;
        else return this->at(key).to_double_or_throw();
      }

      // Attempts to get an integer from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      long get_long(const std::string& key) {
        check_if_object(key);
        if (has_key(key)) return this->at(key).to_long_or_throw();
        else throw marley::Error("Missing JSON key '" + key + '\'');
        return 0.;
      }

      // Attempts to get an integer from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      long get_long(const std::string& key, long default_value) {
        check_if_object(key);
        if (!has_key(key)) return default_value;
        else return this->at(key).to_long_or_throw();
      }

      // Attempts to get a bool from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      bool get_bool(const std::string& key) {
        check_if_object(key);
        if (has_key(key)) return this->at(key).to_bool_or_throw();
        else throw marley::Error("Missing JSON key '" + key + '\'');
        return 0.;
      }

      // Attempts to get a bool from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      bool get_bool(const std::string& key, bool default_value) {
        check_if_object(key);
        if (!has_key(key)) return default_value;
        else return this->at(key).to_bool_or_throw();
      }

      // Attempts to get a string from a JSON object with
      // a given key. If the attempt fails, throw a marley::Error.
      std::string get_string(const std::string& key) {
        check_if_object(key);
        if (has_key(key)) return this->at(key).to_string_or_throw();
        else throw marley::Error("Missing JSON key '" + key + '\'');
        return std::string("");
      }

      // Attempts to get a string from a JSON object with
      // a given key. If the key doesn't exist, use a default value. If a
      // conversion attempt fails, throw a marley::Error.
      std::string get_string(const std::string& key,
        const std::string& default_value)
      {
        check_if_object(key);
        if (!has_key(key)) return default_value;
        else return this->at(key).to_string_or_throw();
      }

      // Copies a subobject from a JSON object with a given key. If the attempt
      // fails, throw a marley::Error, unless the user asks us not to do so.
      marley::JSON get_object(const std::string& key, bool throw_error = true)
      {
        check_if_object(key);
        if (has_key(key)) return this->at(key);
        else if (throw_error) marley::Error("Missing JSON key '" + key + '\'');
        return JSON::make(JSON::DataType::Object);
      }

    private:

      DataType type_ = DataType::Null;
  };

  JSON array() {
    return JSON::make(JSON::DataType::Array);
  }

  template <typename... T> JSON array( T... args ) {
    JSON arr = JSON::make(JSON::DataType::Array);
    arr.append(args...);
    return arr;
  }

  JSON object() {
    return JSON::make(JSON::DataType::Object);
  }

  namespace {

    JSON parse_next(std::istream&);

    void issue_parse_warning(char found_char, const std::string& message)
    {
      std::string warning(message);
      if (found_char == std::ifstream::traits_type::eof())
        warning += "end-of-file";
      else warning += std::string("\'") + found_char + '\'';
      MARLEY_LOG_WARNING() << warning;
    }

    void issue_parse_warning(const std::string& found_str,
      const std::string& message, const std::istream& is)
    {
      std::string warning(message);
      if (!is) warning += "end-of-file";
      else warning += '\'' + found_str + '\'';
      MARLEY_LOG_WARNING() << warning;
    }

    // Skips single-line comments // and multi-line comments /* */
    // These are technically not valid in JSON (the standard doesn't allow
    // comments), but they are valid in Javascript object literals.
    void skip_comment(std::istream& in, bool is_multiline = false) {
      if (is_multiline) {
        char c;
        while (in.get(c)) {
          if (c == '*' && in.peek() == '/') {
            in.ignore();
            break;
          }
        }
      }
      // Ignore all further characters until either a newline or end-of-file
      else in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Skips whitespace and comments, saving the last character read to
    // read_char.
    void skip_ws(std::istream& in, char& read_char) {
      while (read_char = in.get(), std::isspace(read_char)) continue;
      if (read_char == '/') {
        char c = in.peek();
        if (c == '/' || c == '*') {
          read_char = in.get();
          skip_comment(in, c == '*');
          return skip_ws(in, read_char);
        }
      }
    }

    // Removes whitespace and comments from the input stream, putting back
    // the first non-whitespace and non-comment character it finds.
    void consume_ws(std::istream& in) {
      static char next;
      skip_ws(in, next);
      in.putback(next);
    }

    // Removes whitespace and comments from the input stream, returning the
    // first non-whitespace and non-comment character it finds.
    char get_next_char(std::istream& in)
    {
      static char next;
      skip_ws(in, next);
      return next;
    }

    JSON parse_object(std::istream& in) {

      JSON object = JSON::make(JSON::DataType::Object);

      for (;;) {

        consume_ws(in);
        JSON key;

        if ( in.peek() == '}' ) {
          in.ignore();
          return object;
        }
        else if ( in.peek() == '\"' ) {
          key = parse_next(in);
        }
        // The key isn't quoted, so assume it's a single word followed
        // by a colon. Note that vanilla JSON requires all keys to be quoted,
        // but Javascript object literals allow unquoted keys.
        else {
          std::string key_str;
          char c;
          while (in.get(c)) {
            if (c == ':' || std::isspace(c)) {
              in.putback(c);
              break;
            }
            key_str += c;
          }
          key = key_str;
        }

        char next = get_next_char(in);
        if ( next != ':' ) {
          issue_parse_warning(next, "JSON object: Expected colon, found ");
          break;
        }

        consume_ws(in);
        JSON value = parse_next(in);
        object[key.to_string()] = value;

        next = get_next_char(in);
        if ( next == ',' ) continue;
        else if ( next == '}' ) break;
        else {
          issue_parse_warning(next, "JSON object: Expected comma, found ");
          break;
        }
      }

      return object;
    }

    JSON parse_array(std::istream& in) {
      JSON array = JSON::make(JSON::DataType::Array);
      unsigned index = 0;

      for (;;) {

        consume_ws(in);
        if (in.peek() == ']') {
          in.ignore();
          return array;
        }

        array[index++] = parse_next(in);
        consume_ws(in);

        char next = in.get();
        if (next == ',') continue;
        else if (next == ']') break;
        else {
          issue_parse_warning(next, "JSON array: Expected ',' or ']'"
            ", found ");
          return JSON::make(JSON::DataType::Array);
        }
      }

      return array;
    }

    JSON parse_string(std::istream& in) {
      JSON str;
      std::string val;
      for(char c = in.get(); c != '\"' ; c = in.get()) {
        if (c == '\\') {
          switch( in.get() ) {
            case '\"': val += '\"'; break;
            case '\\': val += '\\'; break;
            case '/' : val += '/' ; break;
            case 'b' : val += '\b'; break;
            case 'f' : val += '\f'; break;
            case 'n' : val += '\n'; break;
            case 'r' : val += '\r'; break;
            case 't' : val += '\t'; break;
            case 'u' : {
              val += "\\u" ;
              for(unsigned i = 1; i <= 4; ++i) {
                c = in.get();
                if ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')
                  || (c >= 'A' && c <= 'F')) val += c;
                else {
                  issue_parse_warning(c, "JSON string: Expected hex character"
                    " in unicode escape, found ");
                  return JSON::make(JSON::DataType::String);
                }
              }
              break;
            }
            default: val += '\\'; break;
          }
        }
        else val += c;
      }
      str = val;
      return str;
    }

    //FIXME
    JSON parse_number(std::istream& in, char old) {
      JSON Number;
       std::string val, exp_str;
      char c = old;
      bool isDouble = false;
      long exp = 0;
      for (;;) {
        if ( (c == '-') || (c >= '0' && c <= '9') )
          val += c;
        else if ( c == '.' ) {
          val += c;
          isDouble = true;
        }
        else
          break;
        c = in.get();
      }
      if ( c == 'E' || c == 'e' ) {
        if ( in.peek() == '-' ){ in.ignore(); exp_str += '-';}
        for (;;) {
          c = in.get();
          if ( c >= '0' && c <= '9' )
            exp_str += c;
          else if ( !std::isspace( c ) && c != ',' && c != ']' && c != '}' ) {
            issue_parse_warning(c, "JSON number: Expected a number for"
              " exponent, found ");
            return JSON::make(JSON::DataType::Null);
          }
          else
            break;
        }
        exp = std::stol( exp_str );
      }
      else if ( !std::isspace( c ) && c != ',' && c != ']' && c != '}' ) {
        issue_parse_warning(c, "JSON number: unexpected character ");
        return JSON::make(JSON::DataType::Null);
      }
      in.putback(c);

      if ( isDouble )
        Number = std::stod( val ) * std::pow( 10, exp );
      else {
        if ( !exp_str.empty() )
          Number = std::stol( val ) * std::pow( 10, exp );
        else
          Number = std::stol( val );
      }
      return  Number ;
    }

    JSON parse_bool(std::istream& in, char old) {
      JSON b;
      std::string s(1, old);
      if (old == 't') {
        for (size_t i = 0; i < 3; ++i) s += in.get();
        if (s == "true") b = true;
      }
      else if (old == 'f') {
        for (size_t i = 0; i < 4; ++i) s += in.get();
        if (s == "false") b = false;
      }
      if (b.type() == JSON::DataType::Null) {
        issue_parse_warning("JSON bool: Expected 'true' or 'false', found ",
          s, in);
        return JSON::make(JSON::DataType::Null);
      }
      return b;
    }

    JSON parse_null(std::istream& in) {
      JSON null;
      std::string s(1, 'n');
      for (size_t i = 0; i < 3; ++i) s += in.get();
      if ( s != "null") {
        issue_parse_warning("JSON null: Expected 'null', found ", s, in);
        return JSON::make(JSON::DataType::Null);
      }
      return null;
    }

    JSON parse_next(std::istream& in) {
      char value = get_next_char(in);
      switch(value) {
        case '[' : return parse_array(in);
        case '{' : return parse_object(in);
        case '\"': return parse_string(in);
        case 't' :
        case 'f' : return parse_bool(in, value);
        case 'n' : return parse_null(in);
        default  :
          if ((value <= '9' && value >= '0') || value == '-')
            return parse_number(in, value);
      }
      // Exit gracefully if there was a problem, but complain a bit.
      if (!in) MARLEY_LOG_WARNING() << "Unexpected end of JSON configuration"
        << " file found\n";
      else MARLEY_LOG_WARNING() << "JSON parse: Unknown starting character '"
        << value << "'\n";
      return JSON();
    }
  }

  JSON JSON::load_file(const std::string& filename) {
    std::ifstream in(filename);
    if (in.good()) return load(in);
    else {
      std::cerr << "Couldn't open the file '" << filename << "'\n";
      return JSON::make(JSON::DataType::Null);
    }
  }

  JSON JSON::load(std::istream& in) {
    return parse_next(in);
  }

  JSON JSON::load(const std::string& str) {
    std::stringstream iss(str);
    return load(iss);
  }

}

// Stream operators for JSON input and output using C++ streams
std::ostream& operator<<(std::ostream &os, const marley::JSON& json) {
  os << json.dump_string();
  return os;
}

std::istream& operator>>(std::istream &is, marley::JSON& json) {
  json = marley::JSON::load(is);
  return is;
}
