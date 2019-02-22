r"""A tree implementation for storing sets of p-adic numbers known as a
p-adic tree.

A p-adic tree is a tree that satisfies the following properties:

 - There is a single p-adics used throughout the tree defined by a
   common pAdicBase object stored in each node.

 - Each node is labeled by a tuple of elements of the number field of
   the p-adics. The length of this tuple is the same for each node in
   the tree and is known as the width of the tree. Furthermore the
   elements of this tuple can only be those defined by the method
   :meth:`representatives` of the pAdicBase object defining the
   p-adics.

 - Each node can have nodes that are its children, which will have
   that node as their parent. All children of a node should have
   distinct labels.

 - There exists a unique node that has no parent, called the root of
   the tree.

Furthermore we have the following properties and definitions

 - Each infinite path of the tree starting at the root uniquely
   defines a tuple of p-adic numbers of length equal to the width. If
   the nodes this path passes through, excluding the root, are labeled
   in order by $(a_{1,1}, \ldots, a_{1,n}), (a_{2,1}, \ldots,
   a_{2,n}), ldots$ then the corresponding unique tuple of p-adic
   numbers is ..MATH::

       (a_{1,1} + a_{2,1} \pi + a_{3,1} \pi^2 + \ldots, \ldots,
       a_{1,n} + a_{2,n} \pi + a_{3,n} \pi^2 + \ldots),

   where $\pi$ is the uniformizer of the p-adics associated with this
   tree.

 - We will say that a tuple of p-adic numbers of length the width of
   this tree is in the tree if and only if the corresponding infinite
   branch is in the tree. This makes the tree represent sets of p-adic
   numbers. In this sense all set theoretic operations on p-adics
   trees are defined.

 - The path from the root to a node in the tree defines all tuples of
   p-adic numbers corresponding to infinite paths through that node up
   to a power of the prime of the p-adics. This means we can get a
   tuple of elements of the number field of the p-adics that defines
   this element. We will call this the representative of a node.

 - We will call the length of a path from the root to a node the level
   of that node. Note that the representative of a node of level $n$
   defines all tuples of p-adics numbers corresponding to paths
   through that node up to the prime of the p-adics to the power $n$.

For the implementation of a p-adic tree there are several
classes. First of all there is the class pAdicNode which implements
the basic behavior of a node in a p-adic tree. It makes use of the
class pAdicNodeCollection, which represents a collection of pAdicNodes
that do not all have the same label, to store its children. Using the
class pAdicNodeCollection_inverted which immitates a
pAdicNodeCollection in which all nodes are the root of a tree with all
possible infinite paths but only a few exceptions, one can make trees
out of pAdicNodes that do not represent the empty set.

The second class is pAdicTree, which provides a set like interface to
a p-adic tree defined by pAdicNodes. Furthermore a pAdicTree stores a
list of names, known as variables, of length equal to the with of the
tree. For a tuple of p-adic numbers in this p-adic tree, the i-th
entry will be linked to the i-th variable. Set theoretic operations on
a pAdicTree will keep into account to which value a variable is
linked.

It is important to note that performing operations directly on a
pAdicNode will change the tree it defines. The p-adic tree stored in a
pAdicTree is assumed to be immutable however. Any operations performed
on a pAdicTree will be performed on a copy of the stored p-adic
tree. Even methods that present nodes in the tree will give nodes of a
copy of the original tree. For long computations in which intermediate
steps don't have to be saved it might therefore be useful to work with
the root of the underlying p-adic tree, rather than a pAdicTree
object.

 EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Joey van Langen (2018-07-13): initial version

"""
import weakref

# ****************************************************************************
#       Copyright (C) 2018 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

class pAdicNode(SageObject):
    r"""A node of a p-adic tree.
    
    A p-adic tree is a tree that satisfies the following properties:

     - There is a single p-adics used throughout the tree defined by a
       common pAdicBase object stored in each node.

     - Each node is labeled by a tuple of elements of the number field of
       the p-adics. The length of this tuple is the same for each node in
       the tree and is known as the width of the tree. Furthermore the
       elements of this tuple can only be those defined by the method
       :meth:`representatives` of the pAdicBase object defining the
       p-adics.

     - Each node can have nodes that are its children, which will have
       that node as their parent. All children of a node should have
       distinct labels.

     - There exists a unique node that has no parent, called the root of
       the tree.

    Furthermore we have the following properties and definitions

     - Each infinite path of the tree starting at the root uniquely
       defines a tuple of p-adic numbers of length equal to the width. If
       the nodes this path passes through, excluding the root, are labeled
       in order by $(a_{1,1}, \ldots, a_{1,n}), (a_{2,1}, \ldots,
       a_{2,n}), ldots$ then the corresponding unique tuple of p-adic
       numbers is ..MATH::

           (a_{1,1} + a_{2,1} \pi + a_{3,1} \pi^2 + \ldots, \ldots,
           a_{1,n} + a_{2,n} \pi + a_{3,n} \pi^2 + \ldots),

       where $\pi$ is the uniformizer of the p-adics associated with this
       tree.

     - We will say that a tuple of p-adic numbers of length the width of
       this tree is in the tree if and only if the corresponding infinite
       branch is in the tree. This makes the tree represent sets of p-adic
       numbers. In this sense all set theoretic operations on p-adics
       trees are defined.

     - The path from the root to a node in the tree defines all tuples of
       p-adic numbers corresponding to infinite paths through that node up
       to a power of the prime of the p-adics. This means we can get a
       tuple of elements of the number field of the p-adics that defines
       this element. We will call this the representative of a node.

     - We will call the length of a path from the root to a node the level
       of that node. Note that the representative of a node of level $n$
       defines all tuples of p-adics numbers corresponding to paths
       through that node up to the prime of the p-adics to the power $n$.
    
    EXAMPLE::

        sage: pAdics = pAdicBase(ZZ, 3)
        sage: R = pAdicNode(pAdics=pAdics)
        sage: R.children
        []
        sage: N = pAdicNode(pAdics=pAdics, coefficients=(2,), full=True)
        sage: R.children.add(N)
        sage: R.children
        [p-adic node represented by (2,) with 3 children]
        sage: N.parent() == R
        True
        sage: N.children.list()
        [p-adic node represented by (2,) with 3 children,
         p-adic node represented by (5,) with 3 children,
         p-adic node represented by (8,) with 3 children]

    """
    
    def __init__(self, parent=None, children=None, pAdics=None,
                 coefficients=None, full=False, width=1):
        r"""Initialize this p-adic node.
        
        INPUT:
        
        - ``parent`` -- A pAdicNode or None (default: None) which is
          either the parent of this node or None if it should be the
          root of a p-adic tree.

        - ``children`` -- A pAdicNodeCollection or None (default:
          None). If it is a pAdicNodeCollection it should be the
          collection of children of this node. If it is None it will
          be initialized as an empty or full collection depending on
          the argument ``full``

        - ``pAdics`` -- A pAdicBase object or None (default:None).
          This is the p-adics that should be used for this node.  It
          should be specified unless `parent` is given, in which case
          it can be None. If the `parent` argument is given the
          p-adics will be set to that of the parent node. Note that
          the p-adics should be the same throughout the tree.

        - ``coefficients`` -- A tuple of numbers or None (default:
          None). This is the label of this node.  Note that for
          working correctly the tuple must consist of representatives
          of the residue field of the p-adics of this node, as given
          by the pAdicBase object that defines the p-adics, however
          there is no check to ensure that this is the case.  If this
          argument is None, the coefficients will be set to a tuple of
          zeroes.

        - ``full`` -- A boolean value (default: False). If set to True
          this node will be assumed to be the root of a p-adic tree
          that contains all possible infinite paths. If set to False
          this node will be assumed to have no children at all. If the
          argument `children` is specified this argument will be
          ignored.

        - ``width`` -- A strictly positive integer (default: 1)
          describing the length of the argument `coefficients`. This
          should be the same throughout the tree. It will be ignored
          if either the `coefficients` or the `parent` argument is
          given. If the `parent` argument is given, it will be set to
          the width of the parent. If the `parent` argument is not
          given, but the `coefficients` argument is, it will be set to
          the length of the tuple given as `coefficients`.

        EXAMPLES::

            sage: pAdics = pAdicBase(ZZ, 5)
            sage: N = pAdicNode(pAdics=pAdics); N
            p-adic node represented by (0,) with 0 children
            sage: pAdicNode(parent=N)
            p-adic node represented by (0,) with 0 children
            sage: pAdicNode(pAdics=pAdics, full=True)
            p-adic node represented by (0,) with 5 children
            sage: pAdicNode(pAdics=pAdics, coefficients=(2, 3))
            p-adic node represented by (2, 3) with 0 children
            sage: pAdicNode(pAdics=pAdics, width=3)
            p-adic node represented by (0, 0, 0) with 0 children

        """
        if parent is None:
            self._parent = None
            if pAdics is None or not isinstance(pAdics, pAdicBase):
                raise ValueError("Should specify pAdics if no parent given.")
            self._pAdics = pAdics
            if coefficients is None:
                self.width = width
                self.coefficients = tuple([self.pAdics().zero()]*self.width)
            elif isinstance(coefficients, tuple):
                self.width = len(coefficients)
                self.coefficients = coefficients
            else:
                raise ValueError("The argument coefficients is not a tuple.")
        elif isinstance(parent, pAdicNode):
            self._parent = weakref.ref(parent)
            self._pAdics = parent.pAdics()
            self.width = parent.width
            if coefficients is None:
                self.coefficients = tuple([self.pAdics().zero()]*self.width)
            elif isinstance(coefficients, tuple):
                if len(coefficients) == self.width:
                    self.coefficients = coefficients
                else:
                    raise ValueError("The coefficients argument has the wrong"+
                                     "length.")
            else:
                raise ValueError("The argument coefficients is not a tuple.")
        else:
            raise ValueError("%s is not a valid parent."%(parent,))
            
        if children is None:
            if full:
                children = pAdicNodeCollection_inverted(self)
            else:
                children = pAdicNodeCollection(self)
        if not isinstance(children, pAdicNodeCollection):
            raise ValueError("%s is not a collection of pAdic nodes"%children)
        self.children = children
        self.children.update_parent(self)
        self._update_sub_tree()
        
    def _similar_node(self, other):
        r"""Check whether this and another node are compatible.
                
        Two nodes are compatible iff they both have the same p-adics
        and width.
        
        Note that this function also checks whether the other object
        given is a p-adic node.
        
        INPUT:
        
        - ``other`` -- Any object.
        
        OUTPUT:

        True if `other` is an instance of pAdicNode with the same
        p-adics and width. False in all other cases.

        """
        if not isinstance(other, pAdicNode):
            raise ValueError("%s is not a p-adic node."%other)
        if self.pAdics() != other.pAdics():
            raise ValueError("The p-adics (" + str(self.pAdics()) + " and " +
                             str(other.pAdics()) + ") do not match.")
        if self.width != other.width:
            raise ValueError("The widths of these p-adic nodes do note match.")
            
    def _set_parent(self, parent):
        r"""Set the parent of this node and update subtree accordingly.

        Changes the parent of this node to the given parent if the
        parent is similar according to :func:`_similar_node`.
        It will also reset all values that depend on the hierarchy of
        the nodes for this node and all nodes below it.
        
        .. NOTE::

        Do not use! For internal purposes only.

        INPUT:
        
        - ``parent`` -- A pAdicNode similar to this one.

        """
        self._similar_node(parent)
        self._parent = weakref.ref(parent)
        #reset variables that are not correct after changing the parent
        self._update_sub_tree()
        
    def _update_sub_tree(self):
        r"""Reset hierarchy dependent variables of this node and all below it.
        
        To prevent having to repeat recursive computations, each
        pAdicNode object caches information that has to be calculated
        recursively.  Since this information depends on the structure
        above the node, this information has to be recalculated when
        this structure changes. This method makes sure this happens.

        """
        self._rep = None
        self.children._update_existing_children()
    
    def parent(self):
        r"""Give the parent of this node
        
        OUTPUT:

        A pAdicNode object that is the parent of this node or None if
        it has no parent.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 2)
            sage: N = pAdicNode(pAdics=pAdics, full=True)
            sage: N.parent() == None
            True
            sage: N.children.list()[1].parent() == N
            True

        """
        if self._parent is None:
            return None
        return self._parent()
        
    def is_root(self):
        r"""Determine whether this node is a root.
        
        OUTPUT:

        True if this node is the root of the p-adic tree it is part
        of. False otherwise.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 7)
            sage: N = pAdicNode(pAdics=pAdics, full=True)
            sage: N.is_root()
            True
            sage: N.children.list()[0].is_root()
            False

        """
        return self.parent() is None
        
    def root(self):
        r"""Give the root of the p-adic tree this node is in.
        
        OUTPUT:

        The unique node above this one that has no parent.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 5)
            sage: N = pAdicNode(pAdics=pAdics, full=True)
            sage: N.root() == N
            True
            sage: N.children.list()[0].root() == N
            True

        """
        if self.is_root():
            return self
        else:
            return self.parent().root()
        
    def representative(self):
        r"""A representation of this node as a tuple of numbers.
        
        The representation of this node will be a tuple of zeroes
        if it is the root and it will be
        .. MATH:
        
            r + c * \pi^(l-1)
            
        where `r` is the representation of its parent, `c` is the
        coefficients of this node and `\pi` is the uniformizer
        returned by the p-adics of the tree this node belongs to.
        
        .. NOTE::

        This function is cached.
        
        OUTPUT:

        A tuple of elements of the number field of the p-adics of the
        tree this node is part of, such that all infinite paths from
        the root through this node correspond to tuples of p-adic
        numbers that are equivalent to this tuple modulo $P^n$, where
        $P$ is the prime of the p-adics and $n$ is the level of this
        node.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 7)
            sage: N = pAdicNode(pAdics=pAdics, full=True, width=2)
            sage: N.representative()
            (0, 0)
            sage: N.children.list()[11].representative()
            (4, 1)

        .. SEEALSO::

            :meth:`level`

        """
        if self.is_root():
            return self.coefficients
        if not hasattr(self, "_rep") or self._rep is None:
            self._rep = list(self.parent().representative())
            for i in range(self.width):
                self._rep[i] += (self.coefficients[i] *
                                 self.pAdics().uniformizer()^(self.level()-1))
            self._rep = tuple(self._rep)
        return self._rep
        
    def quotient_tuple(self):
        r"""A representation of this node as a tuple of numbers in a
        corresponding quotient ring.
        
        A representation as a tuple of elements of a quotient ring is
        the same representation as returned by :meth:`representative`
        but then considered as an element of $R / (P^n)$, where $R$ is
        the ring of integers of the p-adics of this tree, $P$ is the
        prime ideal of the p-adics of this tree and $n$ is the level of this
        node.
        
        OUTPUT:

        A tuple of elements of $R / (P^n)$ where $R$ is the ring of
        integers of the p-adics of this node, $P$ is the prime ideal
        of the p-adics of this tree and $n$ is the level of this
        node. This is the reduction modulo $P^n$ of each tuple of
        p-adic numbers that corresonds to an infinite path from the
        root through this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 5)
            sage: N = pAdicNode(pAdics=pAdics, full=True, width=2)
            sage: N.quotient_tuple()
            (0, 0)
            sage: N.quotient_tuple()[0].parent()
            Ring of integers modulo 1
            sage: N2 = N.children.list()[7].children.list()[13]
            sage: N2.quotient_tuple()
            (17, 11)
            sage: N2.quotient_tuple()[0].parent()
            Ring of integers modulo 25

        .. SEEALSO::

            :meth:`representative`
            :meth:`level`

        """
        S = self.pAdics().quotient_ring(self.level())
        return tuple([S(a) for a in self.representative()])
        
    def pAdics(self):
        r"""Return the p-adics of the p-adic tree.
        
        OUTPUT:

        A pAdicBase object that has the p-adics of the p-adic tree
        this node belongs to.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics=pAdics, full=True)
            sage: N.pAdics()
            p-adic information of Rational Field with respect to (3)
            sage: N.children.list()[0].pAdics()
            p-adic information of Rational Field with respect to (3)

        """
        if not hasattr(self, "_pAdics") or self._pAdics is None:
            if self.is_root():
                raise ValueError("A p-adic tree does not have p-adic information.")
            self._pAdics = self.parent().pAdics()
        return self._pAdics
        
    def level(self):
        r"""Give the level of this node.
        
        OUTPUT:

        The length of the shortest path from the root to this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, full=True)
            sage: N = R; N.level()
            0
            sage: N = N.children.list()[0]; N.level()
            1
            sage: N = N.children.list()[0]; N.level()
            2
            sage: N = N.children.list()[0]; N.level()
            3

        """
        if self.is_root():
            return 0
        else:
            return self.parent().level() + 1
            
    def parent_at_level(self, n, my_level=None):
        r"""Give the node at level `n` above this node.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of the node
          we are looking for. This should at most be the level of this
          node.

        - `my_level` -- A non-negative integer (default:
          self.level()).  This argument should equal the level of this
          node.  It is used internally to prevent the recursion needed
          for self.level(). Do not use unless you know the level of
          this node beforehand.
          
        OUTPUT:

        The pAdicNode at level `n` that is in the parent chain of this
        node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, full=True)
            sage: N = R.children_at_level(4)[13]
            sage: N.parent_at_level(0)
            p-adic node represented by (0,) with 2 children
            sage: N.parent_at_level(0) == R
            True
            sage: N.parent_at_level(1)
            p-adic node represented by (1,) with 2 children
            sage: N.parent_at_level(2)
            p-adic node represented by (3,) with 2 children
            sage: N.parent_at_level(3)
            p-adic node represented by (3,) with 2 children
            sage: N.parent_at_level(4)
            p-adic node represented by (11,) with 2 children
            sage: N.parent_at_level(4) == N
            True
            sage: N.parent_at_level(5)
            Traceback (most recent call last):
            ...
            ValueError: p-adic node represented by (11,) with 2 children does not have a parent at level 5

        .. SEEALSO::

            :meth:`level`
            :meth:`parent`
        
        """
        if my_level is None:
            my_level = self.level()
        if n == my_level:
            return self
        elif n < my_level:
            return self.parent().parent_at_level(n, my_level=my_level-1)
        else:
            raise ValueError(str(self) + " does not have a parent at level " +
                             str(n))
    
    def count_children_at_level(self, n, my_level=None):
        r"""Give the number of nodes at level `n` below this one.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of which we
          want to count the number of nodes.

        - `my_level` -- A non-negative integer (default:
          self.level()).  This should always equal the level of this
          node. It is used internally to prevent the recursion needed
          for self.level(). Do not use unless you know the level of
          this node beforehand.
          
        OUTPUT:

        A non-negative integer equal to the number of nodes at level
        `n` that have a path to them from the root through this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 5)
            sage: R = pAdicNode(pAdics=pAdics, full=True)
            sage: N = R.children_at_level(2)[13]
            sage: N.count_children_at_level(4)
            25
            sage: N.count_children_at_level(3)
            5
            sage: N.count_children_at_level(2)
            1
            sage: N.count_children_at_level(1)
            0

        .. SEEALSO::

            :meth:`children_at_level`
            :meth:`level`
            :meth:`parent_at_level`

        """
        if my_level is None:
            my_level = self.level()
        if n < my_level:
            return 0
        elif n == my_level:
            return 1
        elif n == my_level + 1:
            return self.children.size()
        elif self.is_full():
            return self.children.size()^(n - my_level)
        else:
            result = 0
            for child in self.children:
                result += child.count_children_at_level(n, my_level=my_level+1)
            return result
    
    def children_at_level(self, n, my_level=None):
        r"""Give a list of nodes at level `n` below this one.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of the nodes
          we want to get.

        - `my_level` -- A non-negative integer (default:
          self.level()). This should always equal the level of this
          node. It is used internally to prevent the recursion needed
          for self.level(). Do not use unless you know the level of
          this node beforehand.
          
        OUTPUT:

        A list of all nodes at level `n` that have a path to them from
        the root through this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, full=True)
            sage: N = R.children_at_level(2)[3]
            sage: N.children_at_level(4)
            [p-adic node represented by (3,) with 2 children,
             p-adic node represented by (11,) with 2 children,
             p-adic node represented by (7,) with 2 children,
             p-adic node represented by (15,) with 2 children]
            sage: N.children_at_level(3)
            [p-adic node represented by (3,) with 2 children,
             p-adic node represented by (7,) with 2 children]
            sage: N.children_at_level(2)
            [p-adic node represented by (3,) with 2 children]
            sage: N.children_at_level(1)
            []

        """
        if my_level is None:
            my_level = self.level()
        if n == my_level + 1:
            return self.children.list()
        if n > my_level:
            result = []
            for child in self.children:
                 result.extend(child.children_at_level(n, my_level=my_level+1))
            return result
        elif n == my_level:
            return [self]
        else:
            return []
            
    def minimum_full_level(self, my_level=None):
        r"""Give the smallest level at which all nodes below this one are full.

        A node is called full if every possible infinite path from the
        root through that node exist.
        
        INPUT:
        
        - `my_level` -- A non-negative integer (default:
        self.level()). This should always equal the level of this
        node. This is used internally to prevent the recursion needed
        for self.level(). Do not use unless you know the level of this
        node for certain.
        
        OUTPUT:

        The smallest integer $n$, such that all nodes at level $n$ and
        below this node are full. Infinity if no such integer exists.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 3)
            sage: R = pAdicNode(pAdics=pAdics, width=2)
            sage: N1 = pAdicNode(pAdics=pAdics, coefficients=(1, 1))
            sage: R.children.add(N1)
            sage: N2 = pAdicNode(pAdics=pAdics, coefficients=(1, 2), full=True)
            sage: R.children.add(N2)
            sage: R.minimum_full_level()
            2
            sage: N1.minimum_full_level()
            2
            sage: N2.minimum_full_level()
            1

        .. SEE_ALSO::

            :meth:`is_full`,
            :meth:`level`,
            :meth:`children_at_level`

        """
        if my_level is None:
            my_level = self.level()
        if self.is_full():
            return my_level
        if self.is_empty():
            return my_level + 1
        return max([child.minimum_full_level(my_level=my_level+1)
                    for child in self.children])
            
    def is_full(self):
        r"""Determine whether this node is full.

        A node is called full if every possible infinite path from the
        root through that node exist.
        
        OUTPUT:

        True if this node and all nodes below it contain all possible
        children. False otherwise.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 3)
            sage: R = pAdicNode(pAdics=pAdics, width=2)
            sage: N1 = pAdicNode(pAdics=pAdics, coefficients=(1, 1))
            sage: N2 = pAdicNode(pAdics=pAdics, coefficients=(1, 2), full=True)
            sage: R.children.add(N1)
            sage: R.children.add(N2)
            sage: R.is_full()
            False
            sage: N1.is_full()
            False
            sage: N2.is_full()
            True

        .. SEEALSO::

            :meth:`is_empty`

        """
        return self.children.is_full()
        
    def is_empty(self):
        r"""Determine whether this node is empty

        A node is called empty if there is no infinite paths from the
        root through this node.
        
        OUTPUT:

        True if all children of this node are empty. False otherwise.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 3)
            sage: R = pAdicNode(pAdics=pAdics, width=2)
            sage: N1 = pAdicNode(pAdics=pAdics, coefficients=(1, 1))
            sage: N2 = pAdicNode(pAdics=pAdics, coefficients=(1, 2), full=True)
            sage: R.children.add(N1)
            sage: R.children.add(N2)
            sage: R.is_empty()
            False
            sage: N1.is_empty()
            True
            sage: N2.is_empty()
            False

        .. SEE_ALSO::

            :meth:`is_full`

        """
        return self.children.is_empty()
        
    def copy(self):
        r"""Give a copy of this node.
        
        .. NOTE::

        The copy does not have a parent, since the copy process copies
        children. Therefore copying the parent would lead to a double
        copy. A copy of the root would result in a deep copy of this
        tree.
        
        OUTPUT:

        A pAdicNode object that is a copy of this one.

        EXAMPLE::

            sage: pAdics = pAdicBase(ZZ, 5)
            sage: R = pAdicNode(pAdics=pAdics, full=True)
            sage: R.copy()
            p-adic node represented by (0,) with 5 children
            sage: N = R.children.list()[3]
            sage: N.copy()
            p-adic node represented by (3,) with 5 children
            sage: N.copy().parent() == None
            True

        """
        return pAdicNode(pAdics=self.pAdics(),
                         children=self.children.copy(),
                         coefficients=self.coefficients,
                         width=self.width)
                         
    def increase_width(self, n, pAdics=None, coefficients=None):
        r"""Give a new tree with width increased by n.
        
        INPUT:
        
        - ``n`` -- A non negative integer equal to the number of
          additional coefficients in the new width.

        - ``pAdics`` -- A pAdicBase object (default self.pAdics())
          that contains the p-adic information on which the resulting
          p-adic node should be based. This argument should be
          redundant.

        - ``coefficients`` -- A tuple of numbers of length `n`. The
          label of the new node will be the coefficients of this node
          extended with the numbers in this tuple. Note that for
          working correctly the tuple must consist of representatives
          of the residue field of the p-adics of this node, as given
          by the argument `pAdics`, however there is no check to
          ensure that this is the case. If this argument is None, this
          will be initialized as a tuple of zeroes.
          
        OUTPUT:

        A pAdicNode with p-adics as given by the argument `pAdics`,
        coefficients equal to the coefficients of this node extended
        by the given tuple `coefficients`, and as children all
        possible increased width versions of children of this node.

        EXAMPLES::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=2); N
            p-adic node represented by (0, 0) with 0 children
            sage: N.increase_width(1)
            p-adic node represented by (0, 0, 0) with 0 children

        The number of nodes increases accordingly when increasing the
        width::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=2, full=True); N
            p-adic node represented by (0, 0) with 9 children
            sage: N.increase_width(2)
            p-adic node represented by (0, 0, 0, 0) with 81 children

        One can decide what coefficients should fit in the new spots::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=2, full=True); N
            p-adic node represented by (0, 0) with 9 children
            sage: N.increase_width(2, coefficients=(2, 1))
            p-adic node represented by (0, 0, 2, 1) with 81 children

        .. SEEALSO::

            :meth:`decrease_width`
            :meth:`permute_coefficients`

        """
        if pAdics is None:
            pAdics = self.pAdics()
        if coefficients is None:
            coefficients = tuple([pAdics.zero()]*n)
        return pAdicNode(pAdics=pAdics,
                         coefficients=(self.coefficients + coefficients),
                         children=self.children._increase_width(n, pAdics),
                         width=(self.width + n))
                         
    def decrease_width(self, indices, pAdics=None):
        r"""Give a copy of this node with a smaller width.
        
        INPUT:
        
        - ``indices`` -- An iterable object containing distinct
          integers between 0 and the width of this node (0 inclusive,
          width exclusive). These are the indices of the coefficients
          that should be present in the returned node.

        - ``pAdics`` -- A pAdicBase object (default: self.pAdics())
          describing the p-adics of tree that the returned node should
          be part of.
          
        OUTPUT:

        A pAdicNode with the p-adics given by `pAdics` and labeled by
        a tuple of numbers in which the i-th number is the
        indices[i]-th number of the coefficients of this
        node. Furthermore the children of the returned node are those
        obtained from calling :meth:`decrease_width` on the children
        of this node with the same arguments.

        EXAMPLES::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=5); N
            p-adic node represented by (0, 0, 0, 0, 0) with 0 children
            sage: N.decrease_width(range(2))
            p-adic node represented by (0, 0) with 0 children

        The number of children changes accordingly when decreasing the
        width::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=5, full=True); N
            p-adic node represented by (0, 0, 0, 0, 0) with 243 children
            sage: N.decrease_width(range(3))
            p-adic node represented by (0, 0, 0) with 27 children

        We can also choose which indices to keep and in what order::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: N = pAdicNode(pAdics = pAdics, width=5, full=True); N
            p-adic node represented by (0, 0, 0, 0, 0) with 243 children
            sage: N2 = N.children.list()[137]; N2
            p-adic node represented by (2, 0, 0, 2, 1) with 243 children
            sage: N2.decrease_width([1, 4, 0])
            p-adic node represented by (0, 1, 2) with 27 children

        .. SEEALSO::

            :meth:`increase_width`
            :meth:`permute_coefficients`

        """
        if pAdics is None:
            pAdics = self.pAdics()
        coefficients = tuple(self.coefficients[i] for i in indices)
        return pAdicNode(pAdics=pAdics, coefficients=coefficients,
                         children=self.children._decrease_width(indices,
                                                                pAdics),
                         width=len(indices))
                         
    def permute_coefficients(self, permutation, from_root=True):
        """Permute the coefficients of this node.
        
        The permutation will be done in such a way that the i-th entry
        in the new odering will be permutation[i]-th entry of the
        original coefficient.
        
        INPUT:
        
        - ``permutation`` -- a list consisting of the integers 0 to
          the width of this node (0 inclusive, width exclusive), which
          should all appear exactly once. The i-th entry of this list
          should be the index of the coefficient in this node that
          should become the i-th coefficient after permutation.

        - ``from_root`` -- A boolean value (default: True). If set to
          true the coefficients of all nodes in the tree will be
          permuted, whilst the permutation will be limited to this
          node and the ones below it otherwise.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 5)
            sage: R = pAdicNode(pAdics=pAdics, full=True, width=4)
            sage: N = R.children.list()[358]; N
            p-adic node represented by (3, 1, 4, 2) with 625 children
            sage: N.permute_coefficients([1, 3, 0, 2])
            sage: N
            p-adic node represented by (1, 2, 3, 4) with 625 children

        .. SEEALSO::

            :meth:`increase_width`
            :meth:`decrease_width`

        """
        if from_root:
            self.root().permute_coefficients(permutation, from_root=False)
        else:
            self.coefficients = tuple([self.coefficients[permutation[i]]
                                       for i in range(self.width)])
            self._rep = None
            self.children.permute_coefficients(permutation)
    
    def sub_tree(self, children=None):
        r"""Obtain the subtree of this tree containing this node.

        The subtree containing this node is the p-adic tree that
        contains only those nodes in this tree through which there is
        a path from the root that also goes through this node.
        
        INPUT:

        - ``children`` -- A list of nodes that are copies of children
          of this node, or None (default: None). If None the subtree
          will contain all possible children of this node. Otherwise
          the subtree will be limited to only those children given in
          this list.
        
        OUTPUT:

        A pAdicNode that is a copy of the root of the p-adic tree that
        this node is part of. This new tree consists of copies of all
        those nodes in the old tree that have a path through them from
        the root that also passes through this node. If the argument
        `children` was not None, the copy of this node in the new tree
        will have as children precisely the nodes in the list given.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, full=True, width=2)
            sage: N = R.children_at_level(2)[7]; N
            p-adic node represented by (3, 2) with 4 children
            sage: R2 = N.sub_tree()
            sage: R2.children_at_level(2)
            [p-adic node represented by (3, 2) with 4 children]

        """
        if children is None:
            node = self.copy()
        else:
            node = pAdicNode(pAdics=self.pAdics(),
                             coefficients=self.coefficients,
                             width=self.width)
            for child in children:
                node.children.add(child)
        if self.is_root():
            return node
        else:
            return self.parent().sub_tree(children=[node])
        
    def remove(self, child=None):
        r"""Remove a node from the p-adic tree.

        This method removes this node from the tree, or if the
        argument `child` is not None, removes that child of this
        node. If by removing a node a node that is not the root has no
        children left that node will also be removed.
        
        .. NOTE::

        The root of a tree can never be removed and attempting to
        remove it will result in a ValueError.
        
        INPUT:
        
        - ``child`` -- A pAdicNode that is a child of this node, or
          the value None (default: None). If not None that child of
          this node wil be removed. If None this node will be removed.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, full=True, width=2)
            sage: N = R.children_at_level(2)[7]; N
            p-adic node represented by (3, 2) with 4 children
            sage: N.remove()
            sage: R.children_at_level(2)
            [p-adic node represented by (0, 0) with 4 children,
             p-adic node represented by (2, 0) with 4 children,
             p-adic node represented by (0, 2) with 4 children,
             p-adic node represented by (2, 2) with 4 children,
             p-adic node represented by (1, 0) with 4 children,
             p-adic node represented by (3, 0) with 4 children,
             p-adic node represented by (1, 2) with 4 children,
             p-adic node represented by (0, 1) with 4 children,
             p-adic node represented by (2, 1) with 4 children,
             p-adic node represented by (0, 3) with 4 children,
             p-adic node represented by (2, 3) with 4 children,
             p-adic node represented by (1, 1) with 4 children,
             p-adic node represented by (3, 1) with 4 children,
             p-adic node represented by (1, 3) with 4 children,
             p-adic node represented by (3, 3) with 4 children]

        """
        if child is None:
            if self.is_root():
                raise ValueError("Can not remove the root of a tree.")
            else:
                self.parent().remove(child=self)
        else:
            self.children.remove(child)
            if self.is_empty() and not self.is_root():
                self.remove()
    
    def _from_root_list(self):
        r"""Give a path from the root to this node as a list.

        This list consists of all nodes along the path, including the
        root and this node.

        """
        if self.is_root():
            return [self]
        else:
            result = self.parent()._from_root_list()
            result.append(self)
            return result

    def _merge_with_list(self, ls, index=None):
        r"""Merge a list of nodes into this node.

        Part of the implementation of merging. Use the method
        :meth:`merge` to merge nodes.

        INPUT:

        - ``ls`` -- A list of pAdicNodes with the same p-adics and
          width as this node. The list should start with the root of a
          p-adic tree. Each node in this list should be the parent of
          the next node in the list if there is one.

        - ``index`` -- A non-negative integer equal to the level of
          this node plus one. Will be set to that value by
          default. Used to prevent unnecessary recursion in
          meth:`level`

        """
        if index is None:
            index = self.level()+1
        if index >= len(ls):
            self._merge(ls[-1])
        else:
            a_child = None
            if self.children.contains(ls[index].coefficients):
                a_child = self.children.get(ls[index].coefficients)
            else:
                if index >= len(ls)-1:
                    a_child = ls[index].copy()
                    self.children.add(a_child)
                    return None # The child does not have to merge anymore!
                else:
                    a_child = pAdicNode(pAdics=ls[index].pAdics(),
                                        coefficients=ls[index].coefficients,
                                        width=self.width)
                self.children.add(a_child)
            a_child._merge_with_list(ls, index=index+1)
    
    def _merge(self, other):
        r"""Merge another node into this node.
        
        Part of the implementation of merging. Use the method
        :meth:`merge` to merge nodes.

        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        """
        if self.is_full():
            return None # Nothing to add
        for child in other.children:
            if self.children.contains(child.coefficients):
                self.children.get(child.coefficients)._merge(child)
            else:
                self.children.add(child.copy())
        
    def merge(self, other, from_root=True):
        r"""Merge another node into this node.
        
        To merge another node into this node means to add copies of
        all children of that node into this node. If a child with the
        same coefficients already exists in this node the new copied
        child will be merged into that child.

        If the option `from_root` was set to True, will make sure that
        the path from the root to the other node in the other p-adic
        tree exists in this p-adic tree and merge the other node with
        the corresponding node in this tree, rather then this
        node.
        
        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        - ``from_root`` -- A boolean value (default: True). This
          determines whether the two nodes should be merged relative
          to the root or directly.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 3)
            sage: R1 = pAdicNode(pAdics=pAdics, width=2, full=False)
            sage: R2 = pAdicNode(pAdics=pAdics, width=2, full=True)
            sage: N = R2.children_at_level(2)[51]; N
            p-adic node represented by (2, 7) with 9 children
            sage: R1.merge(N)
            sage: R1.children_at_level(2)
            [p-adic node represented by (2, 7) with 9 children]

        .. SEEALSO::

            :meth:`cut`,
            :meth:`limit_to`,
            :meth:`complement`

        """
        self._similar_node(other)
        if from_root:
            self.root()._merge_with_list(other._from_root_list())
        else:
            self._merge(other)
    
    def _cut_list(self, ls, index=None):
        r"""Cut a list of nodes out of this p-adic tree.
        
        Part of the implementation of cutting. Use the method
        :meth:`cut` to cut out nodes.

        INPUT:

        - ``ls`` -- A list of pAdicNodes with the same p-adics and
          width as this node. The list should start with the root of a
          p-adic tree. Each node in this list should be the parent of
          the next node in the list if there is one.

        - ``index`` -- A non-negative integer equal to the level of
          this node plus one. Will be set to that value by
          default. Used to prevent unnecessary recursion in
          meth:`level`

        """
        if index is None:
            index = self.level()+1
        if index >= len(ls):
            self._cut(ls[-1])
        elif self.children.contains(ls[index].coefficients):
            child = self.children.get(ls[index].coefficients)
            child._cut_list(ls, index=index+1)
        # In the remaining case there is nothing to do
    
    def _cut(self, other):
        r"""Cut a node out of this p-adic tree.
        
        Part of the implementation of cutting. Use the method
        :meth:`cut` to cut out nodes.

        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        """
        if (not self.is_root()) and other.is_full():
            self.remove()
        else:
            for child in other.children:
                if self.children.contains(child.coefficients):
                    self.children.get(child.coefficients)._cut(child)
            
    def cut(self, other, from_root=True):
        r"""Cut out a node from this node.
        
        To cut out a node from this node means to remove all infinite
        paths from this node through its children which correspond to
        infinite paths from the other node through its children. In
        practice this means removing all children from this node for
        which their counterparts in the other node are full and
        cutting out those children of the other node that are not full
        from the corresponding children in this node recursively.

        If the argument `from_root` was set to True, will try to find
        the counterpart of the other node in this p-adic tree and cut
        out the other node from that node. By counterpart we mean the
        node in this tree whose representative is the same as the
        representative of the other node. If such a node does not
        exist, this method does nothing.

        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        - ``from_root`` -- A boolean value (default: True). This
          determines whether the other node should be cut out relative
          to the root or directly from this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R = pAdicNode(pAdics=pAdics, width=2, full=True)
            sage: N = R.children_at_level(2)[7]; N
            p-adic node represented by (3, 2) with 4 children
            sage: R.cut(N)
            sage: R.children_at_level(2)
            [p-adic node represented by (0, 0) with 4 children,
             p-adic node represented by (2, 0) with 4 children,
             p-adic node represented by (0, 2) with 4 children,
             p-adic node represented by (2, 2) with 4 children,
             p-adic node represented by (1, 0) with 4 children,
             p-adic node represented by (3, 0) with 4 children,
             p-adic node represented by (1, 2) with 4 children,
             p-adic node represented by (0, 1) with 4 children,
             p-adic node represented by (2, 1) with 4 children,
             p-adic node represented by (0, 3) with 4 children,
             p-adic node represented by (2, 3) with 4 children,
             p-adic node represented by (1, 1) with 4 children,
             p-adic node represented by (3, 1) with 4 children,
             p-adic node represented by (1, 3) with 4 children,
             p-adic node represented by (3, 3) with 4 children]

        .. SEEALSO::

            :meth:`merge`,
            :meth:`limit_to`,
            :meth:`complement`

        """
        self._similar_node(other)
        if from_root:
            self.root()._cut_list(other._from_root_list())
        else:
            self._cut(other)
    
    def _limit_to_list(self, ls, index=None):     
        r"""Limit this p-adic tree to a list of nodes.
        
        Part of the implementation of limiting. Use the method
        :meth:`limit_to` to limit nodes.

        INPUT:

        - ``ls`` -- A list of pAdicNodes with the same p-adics and
          width as this node. The list should start with the root of a
          p-adic tree. Each node in this list should be the parent of
          the next node in the list if there is one.

        - ``index`` -- A non-negative integer equal to the level of
          this node plus one. Will be set to that value by
          default. Used to prevent unnecessary recursion in
          meth:`level`

        """
        if index is None:
            index = self.level() + 1
        if index >= len(ls):
            self._limit_to(self, ls[-1])
        else:
            child = None
            if self.children.contains(ls[index].coefficients):
                child = self.children.get(ls[index].coefficients)
            self.children = pAdicNodeCollection(self)
            if not child is None:
                self.children.add(child)
            
    def _limit_to(self, other):    
        r"""Limit this p-adic tree to a node.
        
        Part of the implementation of limiting. Use the method
        :meth:`limit_to` to limit nodes.

        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        """
        if self.is_full():
            self.children = other.children.copy()
            self.children.update_parent(self)
            return
        removal_list = []
        for child in self.children:
            if other.children.contains(child.coefficients):
                child._limit_to(other.children.get(child.coefficients))
            else:
                removal_list.append(child)
        for child in removal_list:
            child.remove()
    
    def limit_to(self, other, from_root=False):    
        r"""Limit this node to another node.
        
        Limiting this node to another node means that we remove all
        children of this node that are not also children of the other
        node and limit the other children of this node to their
        counterparts in the other node.

        If the argument `from_root` was set to True, will try to find
        the counterpart of the other node in this p-adic tree and
        limit that node to the other node. By counterpart we mean the
        node in this tree whose representative is the same as the
        representative of the other node. If such a node does not
        exist, this method does nothing.
        
        INPUT:

        - ``other`` -- A pAdicNode with the same p-adics and width as
          this node.

        - ``from_root`` -- A boolean value (default: True). This
          determines whether the other node should be cut out relative
          to the root or directly from this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R1 = pAdicNode(pAdics=pAdics, width=2, full=False)
            sage: R2 = pAdicNode(pAdics=pAdics, width=2, full=True)
            sage: N = R2.children_at_level(2)[7]; N
            p-adic node represented by (3, 2) with 4 children
            sage: R1.merge(N)
            sage: R2.limit_to(R1)
            sage: R2.children_at_level(2)
            [p-adic node represented by (3, 2) with 4 children]

        .. SEEALSO::

            :meth:`merge`,
            :meth:`cut`,
            :meth:`complement`

        """
        self._similar_node(other)
        if from_root:
            self.root()._limit_to_list(other._from_root_list())
        else:
            self._limit_to(other)
            
    def complement(self):
        r"""Give the complement of this node.
        
        The complement of a node is a full node with this node cut
        out.

        OUTPUT:

        A pAdicNode that is the complement of this node.

        EXAMPLE::

            sage: pAdics = pAdicBase(QQ, 2)
            sage: R1 = pAdicNode(pAdics=pAdics, width=2, full=True)
            sage: R1.children_at_level(1)[3].remove()
            sage: R1.children_at_level(2)[0].remove()
            sage: R1.children_at_level(2)[4].remove()
            sage: R1.children_at_level(2)[7].remove()
            sage: R1.children_at_level(2)
            [p-adic node represented by (2, 0) with 4 children,
             p-adic node represented by (0, 2) with 4 children,
             p-adic node represented by (2, 2) with 4 children,
             p-adic node represented by (1, 0) with 4 children,
             p-adic node represented by (1, 2) with 4 children,
             p-adic node represented by (3, 2) with 4 children,
             p-adic node represented by (0, 1) with 4 children,
             p-adic node represented by (0, 3) with 4 children,
             p-adic node represented by (2, 3) with 4 children]
            sage: R2 = R1.complement()
            sage: R2.children_at_level(2)
            [p-adic node represented by (0, 0) with 4 children,
             p-adic node represented by (3, 0) with 4 children,
             p-adic node represented by (2, 1) with 4 children,
             p-adic node represented by (1, 1) with 4 children,
             p-adic node represented by (3, 1) with 4 children,
             p-adic node represented by (1, 3) with 4 children,
             p-adic node represented by (3, 3) with 4 children]

        .. SEEALSO::

            :meth:`merge`,
            :meth:`cut`,
            :meth:`limit_to`

        """
        return pAdicNode(pAdics=self.pAdics(), coefficients=self.coefficients,
                         children=self.children.complement(),
                         width=self.width)
        
    def _repr_(self):
        return ("p-adic node represented by " + str(self.representative()) +
                " with " + str(self.children.size()) + " children")
        
    def __eq__(self, other):
        return (self._similar_node(other) and
                self.children == other.children)
               
    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.coefficients, self.children))
        
class pAdicNodeCollection(SageObject):
    r"""
    A collection of p-adic nodes indexed by their coefficients.
    
    A pAdicNodeCollection is a collection of pAdicNode objects
    that never contains two pAdicNode objects with the same
    coefficients. All p-adic nodes in this collection should
    have the same p-adics and width, hence the amount of such
    nodes that can be in this collection is limited.
    
    This class provides a basic implementation of such a
    collection, but similar classes overwrite this class
    and its methods.
    
    For convenience this class can also store the following data
    - A pAdicNode called the parent, which is a node that
      has this collection as its collection of children.
    - A pAdicBase object giving the p-adics shared by the nodes
      in this collection.
    - An strictly positive integer that is the width shared by
      the nodes in this collection.
    """
    
    def __init__(self, parent, pAdics=None, width=1):
        r"""
        Initializes a pAdicNodeCollection object.
        
        INPUT:
        
        - ``parent`` -- A pAdicNode that contains this object as
          its collection of children, or None if no such node
          exists.
        - ``pAdics`` -- A pAdicBase object or None
          (default: None). This gives the p-adics that will be
          shared by the p-adic nodes in this collection. If the
          argument ``parent`` is specified it will revert to the
          p-adics of this parent. Otherwise it is not allowed
          to be None.
        - ``width`` -- A strictly positive integer that is the
          common width among the p-adic nodes in this collection.
        """
        if parent is None and pAdics is None:
            raise ValueError("Should specify pAdics.")
        if isinstance(parent, pAdicNode):
            self._parent = weakref.ref(parent)
            self._pAdics = parent.pAdics()
            self.width = parent.width
        elif parent is None:
            self._parent = None
            if isinstance(pAdics, pAdicBase):
                self._pAdics = pAdics
            else:
                raise ValueError("%s is not a pAdicBase."%(pAdics,))
            if width in ZZ and width >= 0:
                self.width = width
            else:
                raise ValueError("%s is not a valid width."%(width,))
        else:
            raise ValueError("%s is not a pAdicNode"%parent)
        self._dict = {}
        
    def update_parent(self, parent):
        r"""
        Changes the parent of this collection and the nodes therein.
        
        The parent of this collection and of each node in this
        collection will be set to the argument ``parent``.
        
        INPUT:
        
        - ``parent`` -- A pAdicNode that contains this collection
          as its collection of children.
        """
        self._check_pAdic_node(parent)
        self._parent = weakref.ref(parent)
        for v in self._dict.itervalues():
            v._set_parent(parent)
            
    def _update_existing_children(self):
        r"""
        Resets the hierarchy dependant variables of all the nodes in
        this collection and all below it.
        
        .. NOTE:
        For internal purposes only.
        
        To prevent having to repeat recursive computations, each pAdicNode
        object caches information that has to be calculated recursively.
        Since this information depends on the structure above the node,
        this information has to be recalculated when this structure
        changes. This method makes sure this happens.
        """
        for v in self._dict.itervalues():
            v._update_sub_tree()
    
    def _check_pAdic_node(self, node):
        r"""
        Checks whether a given object is a suitable p-adic node.
        """
        if not isinstance(node, pAdicNode):
            raise ValueError("%s is not a pAdicNode"%node)
        if node.pAdics() != self.pAdics():
            raise ValueError("%s does not have p-adics like %s"%(node,
                                                                 self.pAdics()))
        if node.width != self.width:
            raise ValueError("%s does not have width%s"%(node, self.width))
            
    def pAdics(self):
        r"""
        Gives the shared p-adics of the nodes in this collection.
        
        OUTPUT:
        A pAdicBase object that is the p-adics for all nodes in
        this collection.
        """
        return self._pAdics

    def parent(self):
        r"""
        Gives the parent of the nodes in this collection.
        
        OUTPUT:
        A pAdicNode object that is the parent of the nodes
        in this collection or None if they have no parent.
        """
        if self._parent is None:
            return None
        return self._parent()
               
    def add(self, node):
        r"""
        Adds a node to this collection.
        
        If there already exists a node in this collection with
        the same coefficients as the given node, this function
        will raise an error.
        
        INPUT:
        
        - ``node`` -- A pAdicNode with the same p-adics and width
          as shared among nodes of this collection.
        """
        self._check_pAdic_node(node)
        if self._dict.has_key(node.coefficients):
            raise ValueError("A node like %s already exists: %s"%(node,
                                               self._dict[node.coefficients]))
        self._dict[node.coefficients] = node
        if not self.parent() is None:
            node._set_parent(self.parent())
        
    def contains(self, coefficients):
        r"""
        Checks whether this collection contains a node with
        given coefficients.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        
        OUTPUT:
        True - if there exists a node in this collection with
               coefficients equal to the argument ``coefficients``.
        False - otherwise
        """
        return self._dict.has_key(coefficients)
        
    def get(self, coefficients):
        r"""
        Retrieves a node from this collection based on its coefficients.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        
        OUTPUT:
        A pAdicNode object in this collection with coefficients equal
        to the argument ``coefficients``. The function will raise an
        error if no such node exists.
        """
        if not self._dict.has_key(coefficients):
            raise ValueError("No node with coefficients %s exists"&(coefficients,))
        return self._dict[coefficients]
        
    def list(self):
        r"""
        Gives a list of the nodes in this collection.
        
        OUTPUT:
        A list of pAdicNode objects that are in this collection.
        """
        return self._dict.values()
        
    def __iter__(self):
        return self._dict.itervalues()
        
    def size(self):
        r"""
        Gives the amount of nodes in this collection.
        
        OUTPUT:
        A non-negative integer equal to the number of nodes in this
        collection.
        """
        return len(self._dict)
        
    def maximal_size(self):
        r"""
        Gives the maximal amount of nodes in this collection.
        
        OUTPUT:
        A non-negative integer equal to the number of possible
        tuples that can be coefficients for nodes in this
        collection.
        """
        return self.pAdics().size_residue_field()^self.width
        
    def is_full(self):
        r"""
        Determine whether this collection is full.
        
        OUTPUT:
        True - If there is a full node in this collection for
               each possible tuple of coefficients.
        False - Otherwise.
        """
        if self.size() != self.maximal_size():
            return False
        for node in self:
            if not node.is_full():
                return False
        return True
        
    def is_empty(self):
        r"""
        Determine whether this collection is empty.
        
        OUTPUT:
        True - If each node in this collection is empty.
        False - Otherwise.
        """
        if self.size() == 0:
            return True
        for node in self:
            if not node.is_empty():
                return False
        return True
        
    def remove(self, node):
        r"""
        Removes a node from this collection
        
        A node is only removed if it was in the
        collection in the first place.
        
        INPUT:
        
        - ``node`` -- A pAdicNode object that is in
          this collection.
        """
        self._check_pAdic_node(node)
        try:
            if self._dict[node.coefficients] == node:
                self._dict.pop(node.coefficients)
        except KeyError:
            pass
            
    def remove_by_coefficients(self, coefficients):
        r"""
        Removes the node with given coefficients from this
        collection.
        
        A node is only removed if a node with the given
        coefficients was in this collection in the first
        place.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        """
        if coefficients in self._dict:
            self._dict.pop(coefficients)

    def copy(self):
        r"""
        Gives a copy of this collection.
        
        .. NOTE:
        The copy does not have a parent, since the
        copy process copies children. Therefore
        copying the parent would lead to a double
        copy.
        
        OUTPUT:
        A pAdicNodeCollection object that contains
        copies of the nodes in this collection.
        """
        result = pAdicNodeCollection(None, pAdics=self.pAdics(),
                                     width=self.width)
        for n in self:
            result.add(n.copy())
        return result
        
    def complement(self):
        r"""
        Gives the complement of this collection.
        
        OUTPUT:
        A pAdicNodeCollection that contains
        - For each node in this collection its complement.
        - For all possible tuples of coefficients that do not
          have a node in this collection a full node, i.e. a
          node with all possible nodes below it.
        """
        if self.is_full():
            return pAdicNodeCollection(None, pAdics=self.pAdics(),
                                       width=self.width)
        result = pAdicNodeCollection_inverted(None, pAdics=self.pAdics(),
                                              width=self.width)
        for node in self:
            result.get(node.coefficients).cut(node, from_root=False)
        return result
        
    def permute_coefficients(self, permutation):
        r"""
        Permutes the coefficients of nodes.
        
        When called change the coefficients of each node
        in this collection and all below them by making
        the i-th coefficient equal to the the
        permutation[i]-th coefficient of the original
        coefficient tuple. For more information see the
        function :func:`permute_coefficients` of each
        node.
        
        INPUT:
        
        - ``permutation`` -- A list of length equal to
          the common width of nodes in this collection,
          of which each entry is a unique non-negative
          integer smaller than that width. The i-th
          entry of this list should be the index of
          the original entry in the coefficients that
          should become the i-th entry in the new
          coefficients.
        """
        resorted_dict = {}
        for node in self._dict.itervalues():
            node.permute_coefficients(permutation, from_root=False)
            resorted_dict[node.coefficients] = node
        self._dict = resorted_dict
        
    def _increase_width(self, n, pAdics):
        r"""
        Increases the width of all nodes in this collection.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection(None, pAdics=pAdics,
                                     width=self.width+n)
        for cfs in pAdics.representatives(width=n):
            for node in self:
                result.add(node.increase_width(n, pAdics=pAdics,
                                               coefficients=cfs))
        return result
        
    def _decrease_width(self, indices, pAdics):
        r"""
        Limits the coefficients of nodes to certain given indices.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection(None, pAdics=pAdics, width=len(indices))
        for node in self:
            new_node = node.decrease_width(indices, pAdics=pAdics)
            if result.contains(new_node.coefficients):
                result.get(new_node.coefficients).merge(new_node,
                                                          from_root=False)
            else:
                result.add(new_node)
        return result
        
    def _repr_(self):
        result = "["
        firstLoop = True
        for c in self._dict.itervalues():
            if firstLoop:
                firstLoop = False
            else:
                result += ",\n"
            result += str(c)
        result += "]"
        return result
        
    def __eq__(self, other):
        if not isinstance(other, pAdicNodeCollection):
            return False
        if self.pAdics() != other.pAdics():
            return False
        if self.size() != other.size():
            return False
        for c in self._dict:
            if not other.contains(c) or self.get(c) != other.get(c):
                return False
        return True
        
    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(tuple(self._dict.itervalues()))
            
class pAdicNodeCollection_inverted(pAdicNodeCollection):
    r"""
    A pAdicNodeCollection that assumes all possible nodes to exist
    except certain stored exceptions.
    
    A pAdicNodeCollection_inverted is a collection of p-adic
    nodes that assumes all these nodes to be full, i.e.
    containing all possible nodes below them, if nothing is
    specified about them. 
    
    For convenience this class can also store the following data
    - A pAdicNode called the parent, which is a node that
      has this collection as its collection of children.
    - A pAdicBase object giving the p-adics shared by the nodes
      in this collection.
    - An strictly positive integer that is the width shared by
      the nodes in this collection.
    """
    
    def __init__(self, parent, pAdics=None, width=1):
        pAdicNodeCollection.__init__(self, parent, pAdics=pAdics, width=width)
        self._removed = []
    
    def add(self, node):
        r"""
        Adds a node to this collection.
        
        If there already exists a node in this collection with
        the same coefficients as the given node, this function
        will raise an error.
        
        INPUT:
        
        - ``node`` -- A pAdicNode with the same p-adics and width
          as shared among nodes of this collection.
        """
        self._check_pAdic_node(node)
        if node.coefficients not in self._removed:
            raise ValueError("A node like %s already exists."%(node,))
        self._dict[node.coefficients] = node
        self._removed.remove(node.coefficients)
        if not self.parent() is None:
            node._set_parent(self.parent())
        
    def contains(self, coefficients):
        r"""
        Checks whether this collection contains a node with
        given coefficients.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        
        OUTPUT:
        True - if there exists a node in this collection with
               coefficients equal to the argument ``coefficients``.
        False - otherwise
        """
        return coefficients not in self._removed
    
    def get(self, coefficients):
        r"""
        Retrieves a node from this collection based on its coefficients.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        
        OUTPUT:
        A pAdicNode object in this collection with coefficients equal
        to the argument ``coefficients``. The function will raise an
        error if no such node exists.
        """
        if coefficients in self._removed:
            raise ValueError("No node with coefficients %s exists."&(coefficients,))
        if coefficients not in self._dict:
            self._dict[coefficients] = pAdicNode(parent=self.parent(),
                                                 coefficients=coefficients,
                                                 full=True,
                                                 pAdics=self.pAdics(),
                                                 width=self.width)
        return self._dict[coefficients]
        
    def remove(self, node):
        r"""
        Removes a node from this collection
        
        A node is only removed if it was in the
        collection in the first place.
        
        INPUT:
        
        - ``node`` -- A pAdicNode object that is in
          this collection.
        """
        self._check_pAdic_node(node)
        if node.coefficients in self._dict \
                                   and self._dict[node.coefficients] == node:
            self._dict.pop(node.coefficients)
        if node.coefficients not in self._removed:
            self._removed.append(node.coefficients)
            
    def remove_by_coefficients(self, coefficients):
        r"""
        Removes the node with given coefficients from this
        collection.
        
        A node is only removed if a node with the given
        coefficients was in this collection in the first
        place.
        
        INPUT:
        
        - ``coefficients`` -- A tuple of numbers of length equal
          to the common width of nodes in this collection. All
          numbers must be representatives of the residue field
          of the p-adics shared among the nodes in this collection.
        """
        if coefficients in self._dict:
            self._dict.pop(coefficients)
        if coefficients not in self._removed:
            self._removed.append(coefficients)
            
    def _increase_width(self, n, pAdics):
        r"""
        Increases the width of all nodes in this collection.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection_inverted(None, pAdics=pAdics,
                                              width=self.width + n)
        F = pAdics.residue_field()
        M = pAdics.residue_field()**n
        for cfs in pAdics.representatives(width=n):
            for cfs0 in self._removed:
                result.remove_by_coefficients(cfs0 + cfs)
            for node in self._dict.itervalues():
                result._dict[coefficients] = node.increase_width(n, pAdics=pAdics,
                                                                 coefficients=cfs)
        return result
        
    def _decrease_width(self, indices, pAdics):
        r"""
        Limits the coefficients of nodes to certain given indices.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection_inverted(None, pAdics=pAdics,
                                              width=len(indices))
        removal_candidates = []
        for coefficients in self._removed:
            removal_candidates.append(tuple([coefficients[i] for i in indices]))
        s = self.pAdics().size_residue_field()^(self.width - len(indices))
        i = 0
        while i < len(removal_candidates):
            c = removal_candidates.count(removal_candidates[i])
            if c == s:
                result.remove_by_coefficients(removal_candidates[i])
            i += c
        for node in self._dict.itervalues():
            new_node = node.decrease_width(indices, pAdics=pAdics)
            if result.contains(new_node.coefficients):
                result.get(new_node.coefficients).merge(new_node,
                                                          from_root=False)
            else:
                result.add(new_node)
        return result
            
    def list(self):
        r"""
        Gives a list of the nodes in this collection.
        
        OUTPUT:
        A list of pAdicNode objects that are in this collection.
        """
        result = []
        for c in self.pAdics().representatives(width=self.width):
            if c not in self._removed:
                result.append(self.get(c))
        return result
        
    def __iter__(self):
        for c in self.pAdics().representatives(width=self.width):
            if c not in self._removed:
                yield self.get(c)
    
    def size(self):
        r"""
        Gives the amount of nodes in this collection.
        
        OUTPUT:
        A non-negative integer equal to the number of nodes in this
        collection.
        """
        return self.maximal_size() - len(self._removed)
        
    def is_full(self):
        r"""
        Determine whether this collection is full.
        
        OUTPUT:
        True - If there is a full node in this collection for
               each possible tuple of coefficients.
        False - Otherwise.
        """
        if len(self._removed) > 0:
            return False
        for node in self._dict.itervalues():
            if not node.is_full():
                return False
        return True
        
    def is_empty(self):
        r"""
        Determine whether this collection is empty.
        
        OUTPUT:
        True - If each node in this collection is empty.
        False - Otherwise.
        """
        if self.size() == 0:
            return True
        if self.maximal_size() - len(self._removed) - len(self._dict) > 0:
            return False
        for node in self._dict.itervalues():
            if not node.is_empty():
                return False
        return True
        
    def copy(self):
        r"""
        Gives a copy of this collection.
        
        .. NOTE:
        The copy does not have a parent, since the
        copy process copies children. Therefore
        copying the parent would lead to a double
        copy.
        
        OUTPUT:
        A pAdicNodeCollection object that contains
        copies of the nodes in this collection.
        """
        result = pAdicNodeCollection_inverted(None, pAdics=self.pAdics(),
                                              width=self.width)
        for c in self._removed:
            result.remove_by_coefficients(c)
        for (k,v) in self._dict.iteritems():
            result._dict[k] = v.copy()
        return result
        
    def permute_coefficients(self, permutation):
        r"""
        Permutes the coefficients of nodes.
        
        When called change the coefficients of each node
        in this collection and all below them by making
        the i-th coefficient equal to the the
        permutation[i]-th coefficient of the original
        coefficient tuple. For more information see the
        function :func:`permute_coefficients` of each
        node.
        
        INPUT:
        
        - ``permutation`` -- A list of length equal to
          the common width of nodes in this collection,
          of which each entry is a unique non-negative
          integer smaller than that width. The i-th
          entry of this list should be the index of
          the original entry in the coefficients that
          should become the i-th entry in the new
          coefficients.
        """
        pAdicNodeCollection.permute_coefficients(self, permutation)
        resorted_removed = []
        for coefficients in self._removed:
            coefficients = tuple([coefficients[permutation[i]]
                                  for i in range(self.width)])
            resorted_removed.append(coefficients)
        self._removed = resorted_removed
        
    def _repr_(self):
        if self.size() == self.maximal_size():
            return "p-adic node collection excluding no coefficients"
        result = "p-adic node collection excluding coefficients:"
        for c in self._removed:
            result += "\n    " + str(c)
        return result
        
    def __eq__(self, other):
        if not isinstance(other, pAdicNodeCollection):
            return False
        if self.pAdics() != other.pAdics():
            return False
        if self.size() != other.size():
            return False
        for c in self._removed:
            if other.contains(c):
                return False
        for c in self._dict:
            if self._dict[c] != other.get(c):
                return False
        for c in other._dict:
            if self.get(c) != other._dict[c]:
                return False
        return True
        
    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((tuple(self._removed), tuple(self._dict.itervalues())))

class pAdicTree(SageObject):
    r"""
    A tree consisting of p-adic nodes that describe the p-adic
    values of certain variables associated to this tree.
    
    A p-Adic tree consists of the following data.
    - A finite collection of ordered variables.
    - A pAdicBase object that determines the p-adics.
    - A pAdicNode object that is the root of a p-adic tree.
    This node should have the same p-adics as the tree and it
    should have as many coefficients as this tree has variables.
    We consider all nodes below this root to also be part of
    this tree.
    
    If we name the variables of a pAdicTree object
    `x_1, \ldots, x_n`, then the representative of a node of
    level `l` in this tree describes the value of
    `(x_1, \ldots, x_n)` modulo `P^l`, where `P` is the prime
    ideal asssociated to the p-adics of this tree.

    The tree stored in a p-Adic tree should be considered to
    be immutable and this data structure should also be
    treated as such. Any methods in this class that would
    give a different tree perform their operations on a copy
    of this tree.
    
    ..SEEALSO:
        
        :class: pAdicNode
    
    EXAMPLES:
        To be added
    """
    
    def __init__(self, variables=None, prime=None, ring=None, pAdics=None,
                 root=None, width=None, full=True):
        r"""
        Initializes this p-adic tree.
        
        INPUT:
        
        - ``variables`` -- A non-empty iterable object containing
          the names that should be used for the variables of this
          tree. These should be in the order as they should be used.
        - ``prime`` -- A prime ideal (default: pAdics.prime_ideal())
          that defines the p-adics of this p-adic tree. Will be
          ignored if the argument pAdics is defined.
        - ``ring`` -- A ring (default: prime.ring()) that can be
          used to define a pAdicBase object. Will be ignored if
          pAdics is defined.
        - ``pAdics`` -- A pAdicBase object (default: root.pAdics()
          or pAdicBase(ring, prime) if root is not given) that
          defines the p-adics to be used by this tree.
        - ``root`` -- A pAdicNode object
          (default: pAdicNode(pAdics=pAdics, full=full, width=width))
          that is the root of some p-adic tree with the same p-adics
          given before.
        - ``width`` -- A strictly positive integer
          (default: len(variables)) that must equal the width of
          the root.
        - ``full`` -- A boolean value (default: True), that
          determines whether the default value of the argument
          root is a full node (True) or an empty node (False).
        """
        if variables is None:
            raise ValueError("Must specify names of variables")
        if not hasattr(variables, '__iter__'):
            variables = [variables]
        self._variables = tuple(str(v) for v in variables)
        
        if root is None:           
            if pAdics is None:
                if prime is None:
                    raise ValueError("At least the argument prime must be set.")
                prime = Ideal(prime)
                if ring is None:
                    ring = prime.ring()
                pAdics = pAdicBase(ring, prime)
            if width is None:
                width = len(variables)
            root = pAdicNode(pAdics=pAdics, full=full, width=width)
        self._root = root
        
        if width is None:
            width = self._root.width
        if width != self._root.width:
            raise ValueError("The width must be equal to that of the root.")
        if width != len(variables):
            raise ValueError("The width must be equal to the number of variables.")            
                                               
    def _repr_(self):
        return ("pAdic tree over %s\n"%self.pAdics().number_field()) + \
        ("for the prime %s\n"%self.pAdics().prime_ideal()) + \
        ("and the variables %s."%(self._variables,))
    
    def root(self):
        r"""
        Returns the root of this tree.

        NOTE:
        
        To make sure the internal structure of the tree can not be
        modified, this returns a copy of the root rather than the
        root itself.
        
        OUTPUT:
        
        A copy of the unique level 0 pAdicNode in this tree.
        """
        return self._root.copy()
        
    def nodes_at_level(self, level):
        r"""
        Gives a list of the nodes in this tree at a certain level.

        NOTE:

        Will return the nodes as part of a tree that is a copy of
        the one stored in this pAdicTree, to prevent accidental
        modification of the internal structure of this tree.
        
        INPUT:
        
        - `level` -- A non-negative integer that is the level of
          which we want to get the nodes.
          
        OUTPUT:
        
        A tuple consisting of
        - A list of copies of all nodes at level `level` that are
          in this tree.
        - The root of the tree that is a copy of this one. This
          must be assigned a value for the tree to keep existing.
        
        ..SEEALSO:
        
            :func:`pAdicNode.children_at_level`
        """
        T = self.root()
        return T.children_at_level(level), T
        
    def pAdics(self):
        r"""
        Returns the p-adics corresponding to this node.
        
        OUTPUT:
        The pAdicBase object corresponding to this node.
        """
        return self._root.pAdics()
                
    def variables(self):
        r"""
        Gives the variables associated to this p-adic tree.
        
        OUTPUT:
        
        A tuple of names of the variables associated to this
        tree. They are ordered as they are ordered in this tree.
        """
        return self._variables
        
    def add_variables(self, variables):
        r"""
        Adds a (list) of variable(s) to this tree.
        
        Will construct a new tree that is the same
        as this one, but with any new variables given
        added at the end. Any given variable will
        only be added if it did not already exist
        in the current tree. The order of the newly
        added variables will remain the same. Newly
        added variables will be assumed to attain
        any possible value.
        
        INPUT:
        
        - ``variables`` -- An iterable object of
          names for variables. These should be unique
          strings.

        OUTPUT:
  
        A pAdicTree with the following characteristics:
         - The variables of this tree are in this order:
           the variables of this tree in order and then
           any variables given in variables that are not
           already part of this tree in order.
         - For each value associated to the variables
           in this tree there are all the possible values
           for the variables in the returned tree.
        """
        if not hasattr(variables, '__iter__'):
            variables = [variables]
        new_variables = self.variables()
        for var in variables:
            if var not in new_variables:
                new_variables.append(var)
        return pAdicTree(variables=new_variables,
                         root=self.root().increase_width(len(new_variables)))
        
    def remove_variables(self, variables):
        r"""
        Removes a (list of) variable(s) from this tree.
        
        Will return a tree that is this tree with all
        the given variables removed from them. Variables
        are only removed if they were part of this tree
        to begin with. The order of the remaining
        variables will be the same as in this tree.
        
        The nodes in this tree will be changed accordingly,
        assuming that removing these variables corresponds
        to removing the corresponding coefficients in the
        nodes. Multiple nodes with the same coefficients
        will be merged.
        
        ..SEEALSO:
        
            :func:`pAdicNode.decrease_width`
        
        INPUT:
        
        - ``variables`` -- An iterable object containing
          the names of the variables that have to be removed.

        OUTPUT:

        A pAdicTree satisfying the following:
         - Its variables are the same as those in this tree
           except for the ones among the given variables.
           Furthermore their order is the same as it is in
           this tree.
         - The values associated to these variables are all
           those for which there exist a value in this tree
           in which the same variables have that value.
        """
        if not hasattr(variables, '__iter__'):
            variables = [variables]
        all_variables = self.variables()
        variables_to_keep = list(all_variables)
        for var in variables:
            if var in variables_to_keep:
                variables_to_keep.remove(var)     
        indices = [all_variables.index(var) for var in variables_to_keep]
        return pAdicTree(variables=variables_to_keep,
                         root=self.root().decrease_width(indices))
        
    def resort_variables(self, variables):
        r"""
        Resorts the variables in this tree.

        Gives a pAdicTree with the same variables as in this
        tree, except that the ordering is the same as the
        ordering given.

        INPUT:

        - ``variables`` -- A iterable object, containing
          the names of the variables in this tree exactly
          ones and ordered in the preferred ordering.

        OUTPUT:

        A pAdicTree with as variables those given in the
        same order as they are given. These variables are
        given the same values as they had in this tree.
        """
        if not hasattr(variables, '__iter__'):
            variables = [variables]
        variables = list(variables)
        myvars = self.variables()
        if not len(variables) == len(myvars):
            raise ValueError("The given argument variables, does not have the same number of variables as this tree (%d != %d)"%(len(variables), len(myvars)))
        try:
            permutation = [myvars.index(variables[i])
                           for i in range(self._root.width)]
        except ValueError:
            raise ValueError("The variables %s and %s do not correspond 1 to 1."%(variables, myvars))
        return pAdicTree(variables=variables,
                         root=self.root().permute_coefficients(permutation))
        
    def change_variables_to(self, variables, ignore_order=False):
        r"""
        Changes the variables to a given set of variables.

        Performs the actions of adding, removing and reordering
        variables to obtain a pAdicTree with variables matching
        the given argument variables.

        INPUT

        - ``variables`` -- An iterable object of names of
          variables. These are the variables the final tree
          should have in the order it should have them.
        - ``ignore_order`` -- A boolean (default: False)
          indicating whether the order of the variables
          should be ignored.

        OUTPUT:

        A pAdicTree satisfying the following properties:
         - The variables of the returned tree are those
           specified in the argument variables. If
           ignore_order was set to False they are also
           in the same order.
         - The values of the returned tree are such that
           any variable that was not in the original tree
           can have any possible value and any variable
           that was part of the original tree can only have
           values for which there was an appropiate value
           in the original tree.
        """
        if not hasattr(variables, '__iter__'):
            variables = [variables]
        variables = tuple(variables)
        varSet = Set(variables)
        myVarSet = Set(self.variables())
        removeSet = myVarSet - varSet
        result = self
        if removeSet.cardinality() > 0:
            result = result.remove_variables(removeSet)
        addSet = varSet - myVarSet
        if addSet.cardinality() > 0:
            result = result.add_variables(addSet)
        if not ignore_order and self.variables() != variables:
            result = result.resort_variables(variables)
        return result

    def is_empty(self):
        r"""
        Tells whether this tree is empty.

        OUTPUT:

        True - If there exists no infinite branches in
               this tree.
        False - Otherwise
        """
        return self._root.is_empty()

    def is_full(self):
        r"""
        Tells whether this tree is full.

        OUTPUT:

        True - If every possible branch in this tree
               exists.
        False - Otherwise
        """
        return self._root.is_full()

    @cached_method
    def get_values_at_level(self, level):
        r"""
        Gives the values in this tree at a given level.

        INPUT:

        - ``level`` -- A non-negative integer.

        OUTPUT:

        A list which contains for each node at the given level
        its representative.
        """
        result = [node.representative() for node in self._root.children_at_level(level)]
        result.sort()
        return result

    @cached_method
    def give_as_congruence_condition(self):
        m = self._root.minimum_full_level()
        modulus = self.pAdics().prime_ideal()^m
        return self.get_values_at_level(m), modulus
        
    def _check_similar_tree(self, other):
        r"""
        Checks whether two trees have the same pAdics.

        INPUT:

        - ``other`` -- Some object.

        OUTPUT:

        TRUE - If object is a pAdicTree using the same
               pAdicBase as this pAdicTree.
        FALSE - In all other cases.

        """
        if not isinstance(other, pAdicTree):
            raise ValueError("The argument %s is not a pAdicTree."%(other,))
        if self.pAdics() != other.pAdics():
            raise ValueError("The two trees don't have the same pAdics.")
        
    def _get_similar_trees(self, other):
        r"""
        Turns two trees into trees with the same variables.

        INPUT:

        - ``other`` - A pAdicTree using the same pAdicBase
          object as this pAdicTree.

        OUTPUT:

        A tuple of pAdicTree's with the same variables and
        pAdics. The first has the same values for the variables
        as this tree, whilst the second has the same values for
        the variables in the pAdicTree other. The variables
        are the union of the sets of variables of both these
        original trees.
        """
        self._check_similar_tree(other)
        variables = list(self.variables())
        for var in other.variables():
            if var not in variables:
                variables.append(var)
        return self.change_variables_to(variables), \
               other.change_variables_to(variables)
        
    def copy(self):
        r"""
        Makes a copy of this pAdicTree.

        OUTPUT:

        A pAdicTree that contains the same variables
        and values as this pAdicTree.
        """
        return pAdicTree(variables=self.variables(), root=self.root())
        
    def union(self, other):
        r"""
        Gives the union of two pAdicTree's.

        INPUT:

        - ``other`` -- A pAdicTree using the same
          pAdicBase as this pAdicTree.

        OUTPUT:

        A pAdicTree satisfying the following properties
        - The pAdics of that tree are the same as those
          for this tree.
        - The variables of that tree are the union of the
          variables of this tree and the given tree other.
        - For each value in that tree, there is either a
          value in this tree which assigns the same values
          to the same variables or a value in the given tree
          other which assigns the same values to the same
          variables.
        """
        T1, T2 = self._get_similar_trees(other)
        T = T1.root()
        T.merge(T2.root())
        return pAdicTree(variables=T1.variables(), root=T)
        
    def intersection(self, other):
        r"""
        Gives the intersection of two pAdicTree's.

        INPUT:

        - ``other`` -- A pAdicTree using the same
          pAdicBase as this pAdicTree.

        OUTPUT:

        A pAdicTree satisfying the following properties
        - The pAdics of that tree are the same as those
          for this tree.
        - The variables of that tree are the union of the
          variables of this tree and the given tree other.
        - For each value in that tree, there is both a
          value in this tree which assigns the same values
          to the same variables and a value in the given tree
          other which assigns the same values to the same
          variables.
        """
        T1, T2 = self._get_similar_trees(other)
        T = T1.root()
        T.limit_to(T2.root())
        return pAdicTree(variables=T1.variables(), root=T)
        
    def difference(self, other):
        r"""
        Gives the difference of two pAdicTree's.

        INPUT:

        - ``other`` -- A pAdicTree using the same
          pAdicBase as this pAdicTree.

        OUTPUT:

        A pAdicTree satisfying the following properties
        - The pAdics of that tree are the same as those
          for this tree.
        - The variables of that tree are the union of the
          variables of this tree and the given tree other.
        - For each value in that tree, there is  a value in
          this tree which assigns the same values to the
          same variables and not a value in the given tree
          other which assigns the same values to the same
          variables.
        """
        T1, T2 = self._get_similar_trees(other)
        T = T1.root()
        T.cut(T2.root())
        return pAdicTree(variables=T1.variables(), root=T)
        
    def complement(self):
        r"""
        Gives the complement of this pAdicTree.

        OUTPUT:

        A pAdicTree satisfying the following properties
        - The pAdics of that tree are the same as those
          for this tree.
        - The variables of that tree are the same as those
          of this tree.
        - For each value not in this tree, there is a value
          in the returned tree and vice versa.
        """
        return pAdicTree(variables=self.variables(), root=self.root().complement())

    def _cache_key(self):
        return self.variables(), self._root
