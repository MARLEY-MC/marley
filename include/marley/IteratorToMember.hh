#pragma once
#include <iterator>

namespace marley {

  /// @brief Template class that creates an iterator to a class member based on
  /// an iterator to the class object.
  /// @details The IteratorToMember class is based on this
  /// <a href="http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2001/0108/becker/list2.htm">
  /// code</a>.
  /// @tparam It type of the original iterator
  /// @tparam T type pointed to by original iterator
  /// @tparam R type of the member we want to point to
  template<
    typename It,
    typename T,
    typename R
    >
  class IteratorToMember
  {

  public:

    // Some typedefs
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
    IteratorToMember(It from, R T::* memptr)
      : m_it(from), m_memptr(memptr){}

    // Operators *, ->, and [] are first forwarded to the contained
    // iterator, then extract the data member.
    reference operator*() const;
    pointer operator->() const;
    reference operator[](difference_type n) const;

    // All operators that have to do with position are forwarded
    // to the contained iterator.
    IteratorToMember& operator++();
    IteratorToMember operator++(int);
    IteratorToMember& operator--();
    IteratorToMember operator--(int);
    IteratorToMember& operator+=(difference_type n);
    IteratorToMember operator+(difference_type n) const;
    IteratorToMember& operator-=(difference_type n);
    IteratorToMember operator-(difference_type n) const;

    inline difference_type
      operator-(const IteratorToMember<It, T, R>& rhs) const
      { return m_it - rhs.m_it; }

    bool operator==(const IteratorToMember<It, T, R>& rhs) const
      { return m_it == rhs.m_it; }

    bool operator!=(const IteratorToMember<It, T, R>& rhs) const
      { return m_it != rhs.m_it; }

    It m_it;

    protected:
    value_type T::* m_memptr;

  };

  // Member function operators
  template<typename It, typename T, typename R>
    inline typename IteratorToMember<It, T, R>::reference
    IteratorToMember<It, T, R>::operator*() const
    { return (*m_it).*m_memptr; }

  template<typename It, typename T, typename R>
    inline typename IteratorToMember<It, T, R>::pointer
    IteratorToMember<It, T, R>::operator->() const
    { return &((*m_it).*m_memptr); }

  template<typename It, typename T, typename R>
    inline typename IteratorToMember<It, T, R>::reference
    IteratorToMember<It, T, R>::operator[](difference_type n) const
    { return m_it[n].*m_memptr; }

  // Prefix operator++
  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>&
    IteratorToMember<It, T, R>::operator++()
    { ++m_it; return *this; }

  // Prefix operator--
  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>&
    IteratorToMember<It, T, R>::operator--()
    { --m_it; return *this; }

  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>&
    IteratorToMember<It, T, R>::operator+=(difference_type n)
    { m_it += n; return *this; }

  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>&
    IteratorToMember<It, T, R>::operator-=(difference_type n)
    { m_it -= n; return *this; }

  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>
    IteratorToMember<It, T, R>::operator+(difference_type n) const
    { return IteratorToMember<It, T, R>(m_it + n, m_memptr); }

  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>
    IteratorToMember<It, T, R>::operator-(difference_type n) const
    { return IteratorToMember<It, T, R>(m_it - n, m_memptr); }

  /// @brief Postfix operator++ (dummy argument)
  /// @note The standard behavior for this operator is to return the
  /// un-incremented value while incrementing the operand. See discussion
  /// <a href="http://stackoverflow.com/a/3846374/4081973">here</a>.
  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>
    IteratorToMember<It, T, R>::operator++(int)
    {
      auto result = IteratorToMember<It, T, R>(*this);
      ++m_it;
      return result;
    }

  /// @brief Postfix operator-- (see note about the postfix operator++)
  template<typename It, typename T, typename R>
    inline IteratorToMember<It, T, R>
    IteratorToMember<It, T, R>::operator--(int)
    {
      auto result = IteratorToMember<It, T, R>(*this);
      --m_it;
      return result;
    }

  /// @brief Make function for convenient construction.
  template<typename It, typename T, typename R>
    IteratorToMember<It, T, R>
    make_IteratorToMember(It it, R T::* memptr)
    { return IteratorToMember<It, T, R>(it, memptr); }

}
