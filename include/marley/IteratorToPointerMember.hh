/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#pragma once
#include <iterator>

namespace marley {

  /// @brief Template class that creates an iterator to a class member based on
  /// an iterator to a pointer (either bare or smart) to the class object.
  /// @details Ideally this class would inherit from marley::IteratorToMember
  /// and merely override the operator*() and operator->() functions, but such
  /// behavior doesn't appear to be currently possible (6/2016) using template
  /// classes.
  /// @tparam It type of the original iterator
  /// @tparam R type of the member we want to point to
  template< typename It, typename R > class IteratorToPointerMember
  {

  public:

    // Some typedefs
    /// @brief Type of the pointer pointed to by the original iterator
    typedef typename std::iterator_traits<It>::value_type
      T_pointer;
    /// @brief Type of the target of the pointer pointed to by the original
    /// iterator
    /// @details We could have used std::remove_pointer<T_pointer>::type here as well,
    /// but that approach does not work with smart pointers (e.g., std::unique_ptr).
    typedef typename std::remove_reference<
      decltype( *std::declval<T_pointer>() ) >::type T;
    typedef typename std::iterator_traits<It>::iterator_category
      iterator_category;
    typedef typename std::iterator_traits<It>::difference_type
      difference_type;
    typedef R value_type;
    typedef R* pointer;
    typedef R& reference;

    /// @brief Construction from an iterator and a pointer to member.
    /// @param from Iterator to the object
    /// @param memptr Pointer to the member of interest
    IteratorToPointerMember(It from, R T::* memptr) :
      m_it(from), m_memptr(memptr){}

    // Operators *, ->, and [] are first forwarded to the contained
    // iterator, then extract the data member.
    virtual reference operator*() const;
    virtual pointer operator->() const;
    reference operator[](difference_type n) const;

    // All operators that have to do with position are forwarded
    // to the contained iterator.
    IteratorToPointerMember& operator++();
    IteratorToPointerMember operator++(int);
    IteratorToPointerMember& operator--();
    IteratorToPointerMember operator--(int);
    IteratorToPointerMember& operator+=(difference_type n);
    IteratorToPointerMember operator+(difference_type n) const;
    IteratorToPointerMember& operator-=(difference_type n);
    IteratorToPointerMember operator-(difference_type n) const;

    inline difference_type
      operator-(const IteratorToPointerMember<It, R>& rhs) const
      { return m_it - rhs.m_it; }

    bool operator==(const IteratorToPointerMember<It, R>& rhs) const
      { return m_it == rhs.m_it; }
    bool operator!=(const IteratorToPointerMember<It, R>& rhs) const
      { return m_it != rhs.m_it; }

    It m_it;

    protected:
    value_type T::* m_memptr;

  };

  // Member function operators
  template<typename It, typename R>
    inline typename IteratorToPointerMember<It, R>::reference
    IteratorToPointerMember<It, R>::operator*() const
    { return (**m_it).*m_memptr; }

  template<typename It, typename R>
    inline typename IteratorToPointerMember<It, R>::pointer
    IteratorToPointerMember<It, R>::operator->() const
    { return &((**m_it).*m_memptr); }

  template<typename It, typename R>
    inline typename IteratorToPointerMember<It, R>::reference
    IteratorToPointerMember<It, R>::operator[](difference_type n) const
    { return m_it[n].*m_memptr; }

  // Prefix operator++
  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>&
    IteratorToPointerMember<It, R>::operator++()
    { ++m_it; return *this; }

  // Prefix operator--
  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>&
    IteratorToPointerMember<It, R>::operator--()
    { --m_it; return *this; }

  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>&
    IteratorToPointerMember<It, R>::operator+=(difference_type n)
    { m_it += n; return *this; }

  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>&
    IteratorToPointerMember<It, R>::operator-=(difference_type n)
    { m_it -= n; return *this; }

  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>
    IteratorToPointerMember<It, R>::operator+(difference_type n) const
    { return IteratorToPointerMember<It, R>(m_it + n, m_memptr); }

  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>
    IteratorToPointerMember<It, R>::operator-(difference_type n) const
    { return IteratorToPointerMember<It, R>(m_it - n, m_memptr); }

  /// @brief Postfix operator++ (dummy argument)
  /// @note The standard behavior for this operator is to return the
  /// un-incremented value while incrementing the operand. See discussion
  /// <a href="http://stackoverflow.com/a/3846374/4081973">here</a>.
  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>
    IteratorToPointerMember<It, R>::operator++(int)
    {
      auto result = IteratorToPointerMember<It, R>(*this);
      ++m_it;
      return result;
    }

  /// @brief Postfix operator-- (see note about the postfix operator++)
  template<typename It, typename R>
    inline IteratorToPointerMember<It, R>
    IteratorToPointerMember<It, R>::operator--(int)
    {
      auto result = IteratorToPointerMember<It, R>(*this);
      --m_it;
      return result;
    }

  /// @brief Make function for convenient construction.
  template<typename It, typename R>
    IteratorToPointerMember<It, R>
    make_IteratorToPointerMember(It it, R std::iterator_traits<It>
      ::value_type::* memptr)
    { return IteratorToPointerMember<It, R>(it, memptr); }
}
