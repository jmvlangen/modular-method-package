r"""
A tree implementation for storing sets of p-adic numbers

This file contains the classes that are used to store sets
of elements in the completion of a number field at a
finite prime using only elements of the number field and
a tree like structure.

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Joey van Langen (2018-07-13): initial version

"""

# ****************************************************************************
#       Copyright (C) 2018 Joey van Langen <j.m.van.langen@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

def _count_zeroes(iterable):
    r"""
    Counts the number of zeroes in a list.
    
    INPUT:
    
    - ``iterable`` -- A list or other iterable
    item containing numbers.
    
    OUTPUT:
    A non-negative integer equal to the number
    of zeroes in the given list or iterable object.
    """
    return sum(i == 0 for i in iterable)

class pAdicNode(SageObject):
    r"""
    A node of a p-adic tree.
    
    A p-adic node consist of the following data, which
    together constitute a node in a unique p-adic tree:
    - A pAdicBase object that determines the p-adics.
    - A tuple of representatives of the residue field of these p-adics,
      called the coefficients of this node.
    - A reference to the parent of this node if it has one.
    - A pAdicNodeCollection of children of this node.
    
    A p-adic tree consists of p-adic nodes that are linked
    together via parent - child connections. A p-adic tree
    satisfies the following properties:
    - If a node is a child of another node, than that second
      node has the first as its parent.
    - There is a unique node without a parent, called the root.
    - Each node has at most one child with a specific tuple of
      representatives.
    
    Note that each node can recursively find the root of
    the tree through its parents. Furthermore, each node has
    a level, i.e. the number of parents above that node,
    which can again be retrieved recursively.
    
    Using the coefficients p-adic nodes also have an interpretation
    as tuples of numbers. For this let `R`, `P` and `\pi` be respectively
    the maximal order, the prime ideal and a uniformizer of the
    p-adics for this node, then we can represent the node as
    ..MATH::
        
        r + c * \pi^(l-1)
        
    where `r` is a representative of its parent (or zero if it has
    none), `c` is the coefficients of this node, and `l` is the
    level of this node.
    
    Note that all p-adic nodes of a specific level `n` give
    unique representatives of a free module over `R / (P^n)`.
    Also appending such representations up to infinity, will
    result in unique representations of p-adic integers.
    
    EXAMPLES:
    Let us explore the possibilities working with a full tree::
        
        sage: pAdics = pAdicBase(ZZ,3)
        sage: T = pAdicNode(pAdics=pAdics, full=True, width=2)
        sage: for n in T.children_at_level(3):
        ....:     ???
    """
    
    def __init__(self, parent=None, children=None, pAdics=None,
                 coefficients=None, full=False, width=1):
        r"""
        Initializes this p-adic node.
        
        INPUT:
        
        - ``parent`` -- A pAdicNode or None (default: None) which
          is either the parent of this node or None if it should
          be the root of a p-adic tree.
        - ``children`` -- A pAdicNodeCollection or None
          (default: None). If it is a pAdicNodeCollection it
          should be the collection of children of this node. If
          it is None it will be initialized as an empty or full
          connection depending on the argument ``full``
        - ``pAdics`` -- A pAdicBase object or None (default:None).
          This is the p-adics that should be used for this node.
          It should be specified unless ``parent`` is given, in
          which case it can be None. If the ``parent`` argument is given
          the p-adics will be set to that of the parent node.
        - ``coefficients`` -- A tuple of numbers or None
          (default: None). This is the coefficients of this node.
          Note that for working correctly the tuple must consist
          of representatives of the residue field of the p-adics
          of this node, as returned by this node's pAdicBase object,
          however there is no check to ensure that this is the case.
          If this argument is None, the coefficients will be set
          to a tuple of zeroes.
        - ``full`` -- A boolean value (default: False). This
          determines whether this node has all possible nodes
          on a level below it or no children at all, respectively
          True and False. If the argument ``children`` is specified
          this argument will be ignored.
        - ``width`` -- A strictly positive integer (default: 1)
          describing how many numbers the coefficients of this
          node should contain. Note that this argument will be
          ignored if either the ``coefficients`` or the ``parent``
          argument is given. If the ``parent`` argument is given,
          it will revert to the width of the parent. If the
          ``parent`` argument is not given, but the ``coefficients``
          argument is, it will revert to the length of the
          tuple given as ``coefficients``.
          
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
            self._parent = parent
            self._pAdics = parent.pAdics()
            self.width = parent.width
            if coefficients is None:
                self.coefficients = tuple([self.pAdics().zero()]*self.width)
            elif isinstance(coefficients, tuple):
                if len(coefficients) == self.width:
                    self.coefficients = coefficients
                else:
                    raise ValueError("The coefficients argument has the wrong length.")
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
        
    def _check_similar_node(self, other):
        r"""
        Checks whether this and another node are compatible.
        
        ..NOTE:
        For internal purposes only.
        
        Two nodes are compatible if:
        - they have the same width.
        - they have the same p-adics.
        
        Note that this function also checks whether ``other`` is a p-adic node.
        
        INPUT:
        
        - ``other`` -- Any object.
        
        OUTPUT:
        True if ``other`` is an instance of pAdicNode with the same p-adics and
        width. False in all other cases.
        """
        if not isinstance(other, pAdicNode):
            raise ValueError("%s is not a p-adic node."%other)
        if self.pAdics() != other.pAdics():
            raise ValueError("The p-adics (%s and %s) do not match."%(self.pAdics(), other.pAdics()))
        if self.width != other.width:
            raise ValueError("The widths of these p-adic ndoes do note match.")
            
    def _set_parent(self, parent):
        r"""
        Sets the parent of this node and updates subtree accordingly.
        
        .. NOTE:
        Do not use! For internal purposes only.
        
        Changes the parent of this node to the given parent if the parent
        is similar according to :func:`_check_similar_node`.
        It will also reset all values that depend on the hierarchy of the
        nodes for this node ann all nodes below it.

        INPUT:
        
        - ``parent`` -- A pAdicNode similar to this one.
        
        """
        self._check_similar_node(parent)
        self._parent = parent
        #reset variables that are not correct after changing the parent
        self._update_sub_tree()
        
    def _update_sub_tree(self):
        r"""
        Resets hierarchy dependant variables of this node and all below it.
        
        .. NOTE:
        For internal purposes only.
        
        To prevent having to repeat recursive computations, each pAdicNode
        object caches information that has to be calculated recursively.
        Since this information depends on the structure above the node,
        this information has to be recalculated when this structure
        changes. This method makes sure this happens.
        """
        self._rep = None
        self.children._update_existing_children()
    
    def parent(self):
        r"""
        Gives the parent of this node
        
        OUTPUT:
        A pAdicNode object that is the parent of this node
        or None if it has no parent.
        """
        return self._parent
        
    def is_root(self):
        r"""
        Determines whether this node is a root.
        
        OUTPUT:
        True if this node is a root, i.e. it has no parent
        and False otherwise.
        """
        return self.parent() is None
        
    def root(self):
        r"""
        Returns the root of the tree of which node is part.
        
        OUTPUT:
        The unique node above this one that has no parent.
        """
        if self.is_root():
            return self
        else:
            return self.parent().root()
        
    def representative(self):
        r"""
        A representation of this node as a tuple of numbers.
        
        The representation of this node will be a tuple of zeroes
        if it is the root and it will be
        .. MATH:
        
            r + c * \pi^(l-1)
            
        where `r` is the representation of its parent, `c` is the
        coefficients of this node and `\pi` is the uniformizer
        returned by this nodes p-adics.
        
        .. NOTE:
        This function is cached.
        
        OUTPUT:
        A representation as described above.
        """
        if self.is_root():
            return self.coefficients
        if not hasattr(self, "_rep") or self._rep is None:
            self._rep = list(self.parent().representative())
            for i in range(self.width):
                self._rep[i] += self.coefficients[i] \
                              * (self.pAdics().uniformizer()**(self.level()-1))
            self._rep = tuple(self._rep)
        return self._rep
        
    def quotient_tuple(self):
        r"""
        A representation of this node as a tuple of numbers in a
        corresponding quotient ring.
        
        This representation is simply the one returned by
        :func:`representative` modulo `P^l`, where `P` is the
        prime ideal corresponding to the p-adics of this node
        and `l` is the level of this node.
        
        OUTPUT:
        A tuple of elements of `R/(P^l)` where `R`, `P` and `l`
        are respectively the maximal order, prime ideal and
        level corresponding to this node. This tuple is a
        reduction of the one returned by :func:`representative`.
        """
        S = self.pAdics().quotient_ring(self.level())
        return tuple([S(a) for a in self.representative()])
        
    def pAdics(self):
        r"""
        Returns the p-adics corresponding to this node.
        
        OUTPUT:
        The pAdicBase object corresponding to this node.
        """
        if not hasattr(self, "_pAdics") or self._pAdics is None:
            if self.is_root():
                raise ValueError("A p-adic tree does not have p-adic information.")
            self._pAdics = self.parent().pAdics()
        return self._pAdics
        
    def level(self):
        r"""
        Returns the level of this node.
        
        OUTPUT:
        0 if this node is a root, otherwise 1 plus the
        level of its parent.
        """
        if self.is_root():
            return 0
        else:
            return self.parent().level() + 1
            
    def parent_at_level(self, n, my_level=None):
        r"""
        Gives a node above this node that lives at level `n`.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of the
          node we are looking for. This should not be larger than the
          level of this node.
        - `my_level` -- Non-negative integer (default: self.level()).
          This argument should equal the level of this node.
          It is used internally to prevent the recursion needed
          for self.level(). Do not use unless you know the level
          of this node beforehand.
          
        OUTPUT:
        The pAdicNode at level `n` that is in the parent chain of this
        node.
        """
        if my_level is None:
            my_level = self.level()
        if n == my_level:
            return self
        elif n < my_level:
            return self.parent().parent_at_level(n, my_level=my_level-1)
        else:
            raise ValueError("Node %s does not have a parent at level %d"%(self, n))
    
    def count_children_at_level(self, n, my_level=None):
        r"""
        Returns the amount of nodes at level `n` below this one.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of
          which we want to count the amount of nodes. This
          should be at least the level of this node.
        - `my_level` -- A non-negative integer (default: self.level()).
          This should always equal the level of this node. It
          is used internally to prevent the recursion needed for
          self.level(). Do not use unless you know the level
          of this node beforehand.
          
        OUTPUT:
        A non-negative number equal to the number of nodes at level `n`
        that are connected through a child chain to this one.
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
        r"""
        Returns a list of nodes at level `n` below this one.
        
        INPUT:
        
        - `n` -- A non-negative integer that is the level of
          which we want to get the nodes. This should be at
          least the level of this node.
        - `my_level` -- A non-negative integer (default: self.level()).
          This should always equal the level of this node. It
          is used internally to prevent the recursion needed for
          self.level(). Do not use unless you know the level
          of this node beforehand.
          
        OUTPUT:
        A list of all nodes at level `n` that have are connected
        to this node through a child chain.
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
            raise ValueError("Node %s does not have children at level %d"%(self, n))
            
    def minimum_full_level(self, my_level=None):
        r"""
        Returns the smallest level at which each node below
        this one is full.
        
        INPUT:
        
        - `my_level` -- A non-negative inger
        (default: self.level()). This should always
        equal the level of this node. This is used internally
        to prevent the recursion neede for self.level().
        Do not use unless you know the level of this node
        for certain.
        
        OUTPUT:
        A non-negative integer n if all nodes at level n below this one
        are full. Infinity if no node below this one is full, including
        this node itself.
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
        r"""
        Determines whether this node is full.
        
        OUTPUT:
        True if this node and all below it contain all possible
        children. False otherwise.
        """
        return self.children.is_full()
        
    def is_empty(self):
        r"""
        Determines whether this node is empty
        
        OUTPUT:
        True if this node contains no children.
        False otherwise.
        """
        return self.children.is_empty()
        
    def copy(self):
        r"""
        Gives a copy of this node.
        
        .. NOTE:
        The copy does not have a parent, since the
        copy process copies children. Therefore
        copying the parent would lead to a double
        copy. A copy of the root would result in a
        deep copy of this tree.
        
        OUTPUT:
        A pAdicNode object that is a copy of this one.
        """
        return pAdicNode(pAdics=self.pAdics(),
                         children=self.children.copy(),
                         coefficients=self.coefficients,
                         width=self.width)
                         
    def extend_variables(self, n, pAdics=None, coefficients=None):
        r"""
        Returns a new tree with width increased by n.
        
        INPUT:
        
        - ``n`` -- A non negative integer equal to the
          number of additional coefficients you want the
          result to have.
        - ``pAdics`` -- A pAdicBase object (default self.pAdics())
          that contains the p-adic information on which the
          resulting p-adic node should be based. This argument
          should be redundant.
        - ``coefficients`` -- A tuple of numbers of length equal
          to the width of this node plus the given argument ``n``.
          This tuple must consist of representatives of the
          residue field as returned by the pAdicBase object
          given in the previous argument. The first `w` numbers
          should be equal to those stored in this node's
          coefficients, in order. By default this will be set
          to the tuple starting with the coefficients of this
          node extended by zeroes.
          
        OUTPUT:
        A pAdicNode consisting of
        - The p-adics given by the argument ``pAdics``
        - The coefficients given by the argument ``coefficients``.
        - All possible children that can be obtained by extending
          the coefficients of a child of this node.
        """
        if pAdics is None:
            pAdics = self.pAdics()
        return pAdicNode(pAdics=pAdics, coefficients=coefficients,
                         children=self.children._extend_variables(n, pAdics),
                         width=self.width+n)
                         
    def limit_variables(self, indices, pAdics=None):
        r"""
        Gives a copy of this node wherein all coefficients are
        limited to the given indices.
        
        INPUT:
        
        - ``indices`` -- An iterable object containing distinct
          integers between 0 and the width of this node
          (0 inclusive, width exclusive). These are the indices
          of the coefficients that should be present in the
          returned node.
        - ``pAdics`` -- A pAdicBase object (default: self.pAdics())
          describing the p-adic information for the node that
          should be returned. Note that this should always
          be the same as the p-adics of this node.
          
        OUTPUT:
        A pAdicNode consisting of
        - As i-th coefficient the inidices[i]-th coefficient of
          this node.
        - As children all possible children of width len(indices)
          which can be obtained by calling limit_variables
          on a child of this node with the same arguments.
        """
        if pAdics is None:
            pAdics = self.pAdics()
        coefficients = tuple([self.coefficients[i] for i in indices])
        return pAdicNode(pAdics=pAdics, coefficients=coefficients,
                         children=self.children._limit_variables(indices, pAdics),
                         width=len(indices))
                         
    def resort_coefficients(self, permutation, from_root=True):
        """
        Resorts the coefficients of this node.
        
        The resorting will take place such that the i-th coefficient
        of the new ordering will be the permutation[i]-th coefficient
        of the original order.
        
        INPUT:
        
        - ``permutation`` -- a list consisting of the integers 0
          to the width of this node (0 inclusive, width exclusive),
          which should all appear exactly once. The i-th entry of
          this list should be the index of the coefficient in this
          node that should become the i-th coefficient after
          permutation.
        - ``from_root`` -- A boolean value (default: True). If set
          to true the coefficients of all nodes in the tree will be
          resorted, whilst the resorting will be limited to this
          node and the ones below otherwise.
        """
        if from_root:
            self.root().resort_coefficients(permutation, from_root=False)
        else:
            self.coefficients = tuple([self.coefficients[permutation[i]]
                                       for i in range(self.width)])
            self._rep = None
            self.children.resort_coefficients(permutation)
    
    def sub_tree(self, children=None):
        r"""
        Obtains the subtree of this tree that contains this node.
        
        INPUT:
        - ``children`` -- A list of pAdicNode that are children
          of this node, or None (default: None). If None the
          subtree will contain all possible children of this
          node, whilst it will otherwise be only limited to
          the children given in this list.
        
        OUTPUT:
        A pAdicNode that is part of a p-adic tree. This p-adic
        tree consists of copies of nodes that are in either a
        parent or child connection to this node. Besides, if
        the argument ``children`` was specified it will only
        have copies of children of this node specified in that
        list.
        """
        if children is None:
            node = self.copy()
        else:
            node = pAdicNode(pAdics=self.pAdics(),
                             coefficients=self.coefficients,
                             width=self.width)
            for child in self.children:
                node.children.add(child)
        if self.is_root():
            return node
        else:
            return self.parent().sub_tree(children=[node])
        
    def remove(self, child=None):
        r"""
        Removes a node from the tree it is part of.
        
        .. NOTE:
        The root of a tree can never be removed and attempting
        to remove it will result in a ValueError.
        
        INPUT:
        
        - ``child`` -- A pAdicNode that is a child of this node,
          or the value None (default: None). If specified, i.e.
          not None, the specified node will be removed as a child
          from this node. Otherwise, this node will be removed.
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
        r"""
        Gives a list of the nodes between the root and this node (inclusive).
        """
        if self.is_root():
            return [self]
        else:
            result = self.parent()._from_root_list()
            result.append(self)
            return result

    def _merge_with_list(self, ls, index=None):
        r"""
        Merges this node with a list of nodes.
        
        .. NOTE:
        Only for internal use. For merging use func:`merge`
        instead.
        """
        if index is None:
            index = self.level()+1
        if index >= len(ls):
            self._merge(ls[-1])
        else:
            #Needs optimization!
            a_child = None
            if self.children.contains(ls[index].coefficients):
                a_child = self.children.get(ls[index].coefficients)
            else:
                if index >= len(ls)-1:
                    a_child = ls[index].copy()
                    self.children.add(a_child)
                    return None #The child does not have to merge anymore!
                else:
                    a_child = pAdicNode(pAdics=ls[index].pAdics(),
                                        coefficients=ls[index].coefficients,
                                        width=self.width)
                self.children.add(a_child)
            a_child._merge_with_list(ls, index=index+1)
    
    def _merge(self, other):
        r"""
        Merges this node with another node.
        
        .. NOTE:
        Only for internal use. For merging use func:`merge`
        instead.
        """
        if self.is_full():
            return None #Nothing to add
        for child in other.children:
            if self.children.contains(child.coefficients):
                self.children.get(child.coefficients)._merge(child)
            else:
                self.children.add(child.copy())
        
    def merge(self, other, from_root=True):
        r"""
        Merges this node with another node.
        
        Merging a node with this one implies that all nodes
        below the other node that are not present below this
        node, will be added below this node. If merging
        from the root, this will all be done relative to
        the root, meaning that a node will be added at the
        level it existed in the other tree. If not merging
        from the root the node will be added as many levels
        below this node as it was below the other node.
        
        Note that merging from the root might add nodes that
        are not below this one, as it considers everything
        relative to the root and not this node. It will
        however not add any node on a level above that
        of the other node that is not in a parent relationship
        with that node.
        
        INPUT:
        - ``other`` -- A pAdicNode with the same p-adics
          and width as this one.
        - ``from_root`` -- A boolean value (default: True).
          This determines whether the two nodes should be merged
          relative to the root or directly.
        """
        self._check_similar_node(other)
        if from_root:
            self.root()._merge_with_list(other._from_root_list())
        else:
            self._merge(other)
    
    def _cut_list(self, ls, index=None):
        r"""
        Cuts a list of nodes out of the tree of this node.
        
        .. NOTE:
        Only for internal use. For cutting use func:`cut`
        instead.
        """
        if index is None:
            index = self.level()+1
        if index >= len(ls):
            self._cut(ls[-1])
        elif self.children.contains(ls[index].coefficients):
            child = self.children.get(ls[index].coefficients)
            child._cut_list(ls, index=index+1)
    
    def _cut(self, other):
        r"""
        Cuts a node out of the tree of this node.
        
        .. NOTE:
        Only for internal use. For cutting use func:`cut`
        instead.
        """
        if (not self.is_root()) and other.is_full():
            self.remove()
        else:
            for child in other.children:
                if self.children.contains(child.coefficients):
                    self.children.get(child.coefficients)._cut(child)
            
    def cut(self, other, from_root=True):
        r"""
        Cut away a node from the tree of this node.
        
        To cut away a node from the tree of this node,
        we remove all nodes from this tree of which
        their subtrees form a subtree of the subtree
        given by the other node. If this is done
        from the root all subtrees will be taken relative
        to the root. If not from the root all subtrees
        will be taken relative to this and the other node.
        
        Alternatively one could interpret cutting another
        node from this one as removing all (infinite)
        branches from this tree that are (infinite) branches
        in the other tree passing through the other node.
        If doing this from the root all branches originate
        from the root, whilst they originate from this node
        and the other node otherwise.
        
        INPUT:
        - ``other`` -- A pAdicNode with the same p-adics
          and width as this one.
        - ``from_root`` -- A boolean value (default: True).
          This determines whether the node should be cut
          relative to the root or directly from this one.
        """
        self._check_similar_node(other)
        if from_root:
            self.root()._cut_list(other._from_root_list())
        else:
            self._cut(other)
    
    def _limit_to_list(self, ls, index=None):     
        r"""
        Limits the tree of this node to only contain nodes in a list
        
        .. NOTE:
        Only for internal use. For limiting use func:`limit_to`
        instead.
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
        r"""
        Limits the tree of this node to nodes below a given node.
        
        .. NOTE:
        Only for internal use. For limiting use func:`limit_to`
        instead.
        """
        removal_list = []
        for child in self.children:
            if other.children.contains(child.coefficients):
                child._limit_to(other.children.get(child.coefficients))
            else:
                removal_list.append(child)
        for child in removal_list:
            child.remove()
    
    def limit_to(self, other, from_root=False):    
        r"""
        Limits the nodes that can be part of this tree.
        
        Limiting to another node will remove all nodes
        below this node that do not have a counterpart below
        the other node. If doing this from the root the
        counterpart is determined as the node along the
        same path from the root. If doing this not from
        the root is determined by the path relative to
        this respectivelty the other node.
        
        INPUT:
        - ``other`` -- A pAdicNode with the same p-adics
          and width as this one.
        - ``from_root`` -- A boolean value (default: False).
          This determines whether the node should be limited
          relative to the root or directly from this one.
        """
        self._check_similar_node(other)
        if from_root:
            self.root()._limit_to_list(other._from_root_list())
        else:
            self._limit_to(other)
            
    def complement(self):
        r"""
        Gives the complement of this node.
        
        The complement of a node is a node that has all
        nodes below it that are not below the original
        node, i.e. it contains exactly all thoses (infinite)
        branches in a tree that are not part of the original
        tree and pass through this node.
        """
        return pAdicNode(pAdics=self.pAdics(), coefficients=self.coefficients,
                         children=self.children.complement(),
                         width=self.width)
        
    def _repr_(self):
        return "p-adic node represented by " + str(self.representative()) \
                + " with " + str(self.children.size()) + " children."
        
    def __eq__(self, other):
        return isinstance(other, pAdicNode) \
               and self.pAdics() == other.pAdics() \
               and self.children == other.children
               
    def __ne__(self, other):
        return not isinstance(other, pAdicNode) \
               or self.pAdics() != other.pAdics() \
               or self.children != other.children
        
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
            self._pAdics = parent.pAdics()
            self.width = parent.width
        elif parent is None:
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
        self.parent = parent
        
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
        self.parent = parent
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
        if not self.parent is None:
            node._set_parent(self.parent)
        
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
        
    def resort_coefficients(self, permutation):
        r"""
        Permutes the coefficients of nodes.
        
        When called change the coefficients of each node
        in this collection and all below them by making
        the i-th coefficient equal to the the
        permutation[i]-th coefficient of the original
        coefficient tuple. For more information see the
        function :func:`resort_coefficients` of each
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
            node.resort_coefficients(permutation, from_root=False)
            resorted_dict[node.coefficients] = node
        self._dict = resorted_dict
        
    def _extend_variables(self, n, pAdics):
        r"""
        Increases the width of all nodes in this collection.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection(None, pAdics=pAdics,
                                     width=self.width+n)
        F = pAdics.residue_field()
        for node in self:
            for m in pAdics.representatives(width=n):
                coefficients = tuple(list(node.coefficients) +
                                     [c for c in m])
                result.add(node.extend_variables(n, pAdics=pAdics,
                                                 coefficients=coefficients))
        return result
        
    def _limit_variables(self, indices, pAdics):
        r"""
        Limits the coefficients of nodes to certain given indices.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection(None, pAdics=pAdics, width=len(indices))
        for node in self:
            new_node = node.limit_variables(indices, pAdics=pAdics)
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
        if not self.parent is None:
            node._set_parent(self.parent)
        
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
            self._dict[coefficients] = pAdicNode(parent=self.parent,
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
            
    def _extend_variables(self, n, pAdics):
        r"""
        Increases the width of all nodes in this collection.
        
        .. NOTE:
        For internal purposes only.
        """
        result = pAdicNodeCollection_inverted(None, pAdics=pAdics,
                                              width=self.width + n)
        F = pAdics.residue_field()
        M = pAdics.residue_field()**n
        for cfs in self._removed:
            for m in M:
                coefficients = tuple(list(cfs) + [F.lift(c) for c in m])
                result.remove_by_coefficients(coefficients)
        for node in self._dict.itervalues():
            for m in M:
                coefficients = tuple(list(node.coefficients) +
                                     [F.lift(c) for c in m])
                result._dict[coefficients] = node.extend_variables(n, pAdics=pAdics,
                                                         coefficients=coefficients)
        return result
        
    def _limit_variables(self, indices, pAdics):
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
            new_node = node.limit_variables(indices, pAdics=pAdics)
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
        
    def resort_coefficients(self, permutation):
        r"""
        Permutes the coefficients of nodes.
        
        When called change the coefficients of each node
        in this collection and all below them by making
        the i-th coefficient equal to the the
        permutation[i]-th coefficient of the original
        coefficient tuple. For more information see the
        function :func:`resort_coefficients` of each
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
        pAdicNodeCollection.resort_coefficients(self, permutation)
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
        self._variables = tuple(variables)
        
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
        
        A list of all nodes at level `level` that are in this tree.
        
        ..SEEALSO:
        
            :func:`pAdicNode.children_at_level`
        """
        return self.root().children_at_level(level)
        
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
                         root=self.root().extend_variables(len(new_variables)))
        
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
        
            :func:`pAdicNode.limit_variables`
        
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
                         root=self.root().limit_variables(indices))
        
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
                         root=self.root().resort_coefficients(permutation))
        
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

    @cached_method
    def give_as_congruence_condition(self):
        m = self._root.minimum_full_level()
        values = [node.representative() for node in self._root.children_at_level(m)]
        modulus = self.pAdics().prime_ideal()^m
        values.sort()
        return values, modulus
        
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
        self._check_similar_tree(self, other)
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
        return pAdicTree(variables=T1.variables(), root=T1.root().merge(T2.root()))
        
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
        return pAdicTree(variables=T1.variables(), root=T1.root().limit_to(T2.root()))
        
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
        return pAdicTree(variables=T1.variables(), root=T1.root().cut(T2.root()))
        
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
