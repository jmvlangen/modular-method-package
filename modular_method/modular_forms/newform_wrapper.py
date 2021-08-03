r"""Tools to work with newforms from different sources.

This file provides a wrapper class that has subclasses wrapping around
different kinds of newforms to provide a uniform way of working with
newforms from different sources.

The base class for wrapped newforms is :class:`WrappedNewform`. It
has subclasses wrapping around a newform produced by Sage
(:class:`WrappedNewform_sage`), wrapping around a newform produced by
Magma (:class:`WrappedNewform_magma`) and wrapping around a newform
produced by reading fourier coefficients from a file
(:class:`WrappedNewform_stored`).

This file also contains methods for saving and loading wrapped
newforms. To store a newform we store its level, the corresponding
character and some of its fourier coefficients. Saving can be done
with the function :func:`save_newforms` and loading can be done with
the function :func:`load_newforms`. The files storing newforms are
human readable.

This file also gives a single method to compute newforms and returned
them as wrapped version, called :func:`get_newforms`. This method has
the option of choosing between Sage, magma and loading from a file to
compute the newforms.

EXAMPLES:

Working with newforms from Sage::

    sage: from modular_method.modular_forms.newform_wrapper import get_newforms
    sage: eps = DirichletGroup(16).gens()[1]
    sage: nf = get_newforms(16, character=eps)[0]; nf
    q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)
    sage: nf.level()
    16
    sage: nf.character()
    Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

Working with newforms from magma::

    sage: from modular_method.modular_forms.newform_wrapper import get_newforms
    sage: eps = DirichletGroup(16).gens()[1]
    sage: nf = get_newforms(16, character=eps, algorithm='magma')[0]; nf
    q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + O(q^12)
    sage: nf.level()
    16
    sage: nf.character()
    Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

Saving and reloading newforms::

    sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
    sage: eps = DirichletGroup(16).gens()[1]
    sage: f = tmp_filename(ext='.nfs')
    sage: save_newforms(get_newforms(16, character=eps), f)
    sage: nf = load_newforms(f)[0]; nf
    q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + (-2*a - 2)*q^12 + (a - 1)*q^13 + (2*a - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-a + 1)*q^18 + (-3*a + 3)*q^19 + O(q^20)
    sage: nf.level()
    16
    sage: nf.character()
    Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

Computations that can be done with newforms::

    sage: from modular_method.modular_forms.newform_wrapper import get_newforms
    sage: nf = get_newforms(19)[0]
    sage: nf.trace_of_frobenius(5)
    3
    sage: nf.determinant_of_frobenius(11, power=2)
    121
    sage: nf.characteristic_polynomial(7)
    x^2 + x + 7

AUTHORS:

- Joey van Langen (2019-03-04): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import re

from sage.structure.sage_object import SageObject

from sage.all import Integer, QQ, ZZ, Integers

from sage.misc.cachefunc import cached_method

from sage.interfaces.magma import magma

from sage.rings.number_field.number_field import is_NumberField, NumberField
from sage.rings.ideal import is_Ideal
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.modular.dirichlet import is_DirichletCharacter, DirichletGroup

from sage.modular.modform.constructor import Newforms

from modular_method.number_fields.dirichlet_characters import dirichlet_product
from modular_method.number_fields.field_constructors import common_embedding_field
from modular_method.number_fields.field_constructors import composite_field

from modular_method.polynomial.symmetric_polynomials import polynomial_to_symmetric

def get_newforms(level, character=None, algorithm='sage', base_field=QQ,
                 names='a', path=None):
    r"""Compute the newforms of weight 2 of a given level and character.

    INPUT:

    - ``level`` -- A non-zero ideal of the ring of integers of the
      number field given as ``base_field`` which is the level of the
      newforms. If ``base_field`` is the rationals then this may also
      be given as a strictly positive integer.

    - ``character`` -- A dirichlet character of which the conductor
      divides the given level.

    - ``algorithm`` -- One of the following possible arguments: 'sage'
      (default) if sage should be used to compute the newforms;
      'magma' if magma should be used to compute the newforms; or
      'file' if the newforms should be loaded from a file.

    - ``base_field`` -- A number field or the rationals (default: QQ)
      over which we want to compute newforms. This determines what
      kind of modular forms are computed. The possibilities are:

       * classical modular forms for base field the rationals
       * Hilbert modular forms for totally real base fields (not
         supported by a sage algorithm)
       * Bianchi modular forms for complex quadratic fields (not
         supported by a sage algorithm)

    - ``names`` -- An argument required by the sage implementation of
      newforms to be used as the names for the generator of
      coefficient fields of newforms that are not QQ.

    - ``path`` -- A string or None (default: None). Only used in case
      the algorithm is set to file, in which case it should be the
      path to the file from which to load the newforms as a string.

    OUTPUT:

    A list of instances of WrappedNewform that contains exactly one
    newform in each galois orbit of newforms of level `level`,
    (parallel) weight 2, base field `base_field`, and character
    `character`.

    EXAMPLES:

    Getting newforms using Sage's algorithm::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: get_newforms(26, algorithm='sage')
        [q - q^2 + q^3 + q^4 - 3*q^5 + O(q^6), q + q^2 - 3*q^3 + q^4 - q^5 + O(q^6)]
        sage: eps = DirichletGroup(16).gens()[1]
        sage: get_newforms(16, character=eps, algorithm='sage')
        [q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)]
        sage: get_newforms(23, algorithm='sage', names='b')
        [q + b0*q^2 + (-2*b0 - 1)*q^3 + (-b0 - 1)*q^4 + 2*b0*q^5 + O(q^6)]

    Using magma to do the computations instead::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: get_newforms(26, algorithm='magma') # optional - magma
        [q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + O(q^12),
         q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 + O(q^12)]
        sage: eps = DirichletGroup(16).gens()[1]
        sage: get_newforms(16, character=eps, algorithm='magma') # optional - magma
        [q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + O(q^12)]

    Getting newforms from a file::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(get_newforms(26), f)
        sage: get_newforms(26, algorithm='file', path=f)
        [q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + q^12 + q^13 + q^14 - 3*q^15 + q^16 - 3*q^17 + 2*q^18 + 2*q^19 + O(q^20),
         q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 - 3*q^12 - q^13 + q^14 + 3*q^15 + q^16 - 3*q^17 + 6*q^18 + 6*q^19 + O(q^20)]

    """
    if algorithm == 'sage':
        if base_field != QQ:
            raise NotImplementedError("Sage algorithm not implemented for " +
                                      "base field " + str(base_field))
        if character is None:
            nfs = Newforms(level, names=names)
        else:
            eps = character.primitive_character().extend(level)
            nfs = Newforms(eps, names=names)
        result = [WrappedNewform_sage(f) for f in nfs]
    elif algorithm == 'magma':
        if base_field == QQ:
            if character is None:
                cfs = magma.CuspForms(level)
                nfs = magma.Newforms(cfs)
            else:
                eps = character.primitive_character()
                gens = eps.parent().unit_gens()
                Dm = magma.DirichletGroup(eps.modulus(), magma(eps.base_ring()))
                candidate = False
                for eps_m in Dm.Elements():
                    candidate = all(eps_m(n) == eps(n) for n in gens)
                    if candidate:
                        break
                if candidate: # We found a matching character
                    eps_m = magma.DirichletGroup(level,
                                                 magma(eps.base_ring()))(eps_m)
                    cfs = magma.CuspForms(eps_m)
                    nfs = magma.Newforms(cfs)
                else:
                    raise ValueError("There is no dirichlet character in magma " +
                                     "matching %s"%(eps,))
            result = [WrappedNewform_magma(orbit[Integer(1) ]) for orbit in nfs]
        elif base_field.is_totally_real():
            if character != None:
                raise NotImplementedError("Magma algorithm not implemented for " +
                                          "Hilbert modular forms with character")
            K = magma(base_field)
            N = magma(level)
            if not N.Parent().Ring().IsMaximal().sage():
                raise ValueError("The given ideal is not of the ring of integers")
            cuspspace = K.HilbertCuspForms(N)
            newspace = cuspspace.NewSubspace()
            newdecomp = newspace.NewformDecomposition()
            result = [WrappedNewform_magma_hilbert(newspace) for newspace in newdecomp]
        elif base_field.degree() == Integer(2) :
            if character != None:
                raise NotImplementedError("Magma algorithm not implemented for " +
                                          "Bianchi modular forms with character")
            K = magma(base_field)
            N = magma(level)
            if not N.Parent().Ring().IsMaximal().sage():
                raise ValueError("The given ideal is not of the ring of integers")
            cuspspace = K.BianchiCuspForms(N)
            newspace = cuspspace.NewSubspace()
            newdecomp = newspace.NewformDecomposition()
            result = [WrappedNewform_magma_bianchi(newspace) for newspace in newdecomp]
        else:
            raise NotImplementedError("Magma algorithm not implemented for " +
                                      "base field " + str(base_field))
    elif algorithm == 'file':
        if character is None:
            character = DirichletGroup(Integer(1) )[Integer(0) ]
        character = character.primitive_character()
        if path is None:
            raise ValueError("Argument path should be set if algorithm file " +
                             "is chosen.")
        to_do = [load_newforms(path)]
        result = []
        while len(to_do) > Integer(0) :
            check = to_do
            to_do = []
            for element in check:
                if isinstance(element, list):
                    to_do.extend(element)
                elif (isinstance(element, WrappedNewform) and
                      ((element.base_field() == QQ and base_field == QQ) or
                      (element.base_field().absolute_degree() == base_field.absolute_degree() and
                       element.base_field().absolute_polynomial() == base_field.absolute_polynomial()))
                      and element.level() == level and
                      element.character().primitive_character() == character):
                    result.append(element)
                else:
                    pass # Not a right newform, skip!
    else:
        raise ValueError("%s is not a valid algorithm."%(algorithm,))
    return result

class LabeledElement:
    r"""A helper class that stores an element and a label."""
    def __init__(self, label, element):
        self.label = label
        self.element = element

def save_newforms(newforms, file_name, coefficients=Integer(50) , repr_coefficients=True,
                  save_cm=True):
    r"""Save newforms to a file.

    Save a newform or a list of newforms to a file. This file will
    store for each newform information about its level, its character
    and some fourier coefficients. It can also store whether or not
    the newform has complex multiplication, but this is optional.

    The file is written in such a format that it should be human
    readable, using whitespace to lay out the file in a more readable
    way. Independent of the whitespace the file can be read again by
    the function :func:`load_newforms` to get again a list of wrapped
    newforms.

    In the file in which the newforms are save we use the following
    notation, written here as regular expressions:

    <list> := '[' ( <element> ( ',' <element> )* )? ']'
    <element> := ( '<' <identifier> '>' ':=' )? ( <list> | <rational> )
    <identifier> := <letter>+
    <rational> := <integer> ( '/' <positive_integer> )?
    <integer> := ( '-' )? <zero> | <positive_integer>
    <positive_integer> := <non_zero_digit> ( <digit> )*
    <digit> := <zero> | <non_zero_digit>
    <non_zero_digit> := [1-9]
    <zero> := '0'
    <letter> := [a-zA-Z]

    Note that for <list>, <element> and <identifier> whitespace
    between the different building blocks is ignored. Furthermore we
    have the following ways of representing different bits of data.

    - A boolean value is represented by the integer 1 if it is True
      and the integer 0 if it is False.

    - An element of a number field is represented as the list of
      rational coefficients with respect to the power basis in the
      generator, preceded by the identifier 'element'

    - A polynomial with rational coefficients is represented by a list
      of its coefficients (starting at the constant term) preceded by
      the identifier 'polynomial'

    - A number field is represented by a list containing a polynomial,
      preceded by the identifier 'field'. The polynomial is the
      defining polynomial of the number field.

    - A (fractional) ideal of a number field is represented by a list
      containing the elements of the number field that generate the
      ideal, preceded by the identifier 'ideal'

    - A list of values of a function is represented by a list
      containing lists of exactly two elements, all preceded by the
      identifier 'values'. The function maps the first element of a
      list in the corresponding list to the second element thereof.

    - A character is represented by a list containing an integer
      preceded by the identifier 'conductor' and a list of values of
      the character preceded by the identifier 'values', all preceded
      by the identifier 'character'. The integer with identifier
      'conductor' will be the conductor of the character. The entry
      labeled values will be pairs of integers, such that if $\zeta$
      is the relevant $n$-th root of unity a pair $(k, e)$ appears in
      this list if the character takes the value $zeta^e$ at $k$.

    - A newform is represented by a list containing an integer or an
      ideal preceded by the identifier 'level', a boolean preceded by
      the identifier 'cm', a character, a number field preceded by the
      identifier 'basefield', a number field preceded by the
      identifier 'coefficientfield' and a list of values preceded by
      the identifier 'value', all preceded by the identifier
      'newform'. The first entry will be the level of the newform, the
      second a boolean indicating whether or not this newform has
      complex multiplication, the third the corresponding character,
      the fourth the base field of the newform, the fifth the
      coefficient field of the newform and the last the coefficients
      of the newform at (some) indices and the traces of frobenius at
      some prime ideals. For backwards compatibility the label
      basefield may be left out (in which case it is interpreted to be
      the rationals) and the label 'coefficientfield' may be replaced
      with an entry that is just a number field. The entry with label
      'cm' may be left out or set to -1 to indicate that this
      information is not known.

    - A list of things is represented as a list of the corresponding
      representations.

    INPUT:

    - ``newforms`` -- An instance of WrappedNewform. This may also be
      a list or other iterable containing as elements instances of
      WrappedNewform or lists that satisfy the same property. These
      are the newforms that will be saved to the file.

    - ``file_name`` -- A string containing the file name to which the
      given newforms should be saved.

    - ``coefficients`` -- An iterable object of non-negative integers
      or a non-negative integer (default: 50). This determines the
      coefficients and traces of frobenius of the newform that will be
      saved. If it is an iterable object and the newform is a
      classical newform, will save all the coefficients with the
      indices given in that object. If it is a single non-negative
      integer will make this to be the list of all non-negative
      integers up to the given integer (excluding the integer
      itself). If the newform is not a classical newform, will instead
      save the traces of frobenius of every prime ideal above a prime
      number in this list.

    - ``repr_coefficients`` -- A boolean value (default: True).  If
      set to true, will always save the first 20 coeficients of each
      classical newform. This is recommmended to have a string
      representation of the newform as a q-expansion.

    - ``save_cm`` -- A boolean value (default: True) indicating
      whether for each newform saved the information whether or not it
      has complex multiplication should also be computed and saved. If
      set to False or whether or not the newform has CM can not be
      determined, this will not be done and the field 'cm' of a
      newform will be set to -1.

    EXAMPLES::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(26); nfs
        [q - q^2 + q^3 + q^4 - 3*q^5 + O(q^6), q + q^2 - 3*q^3 + q^4 - q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)
        [q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + q^12 + q^13 + q^14 - 3*q^15 + q^16 - 3*q^17 + 2*q^18 + 2*q^19 + O(q^20),
         q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 - 3*q^12 - q^13 + q^14 + 3*q^15 + q^16 - 3*q^17 + 6*q^18 + 6*q^19 + O(q^20)]

    One can store multiple newform lists in a single file::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs1 = get_newforms(26); nfs1
        [q - q^2 + q^3 + q^4 - 3*q^5 + O(q^6), q + q^2 - 3*q^3 + q^4 - q^5 + O(q^6)]
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nfs2 = get_newforms(16, character=eps); nfs2
        [q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms([nfs1, nfs2], f)
        sage: load_newforms(f)
        [[q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + q^12 + q^13 + q^14 - 3*q^15 + q^16 - 3*q^17 + 2*q^18 + 2*q^19 + O(q^20),
          q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 - 3*q^12 - q^13 + q^14 + 3*q^15 + q^16 - 3*q^17 + 6*q^18 + 6*q^19 + O(q^20)],
         [q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + (-2*a - 2)*q^12 + (a - 1)*q^13 + (2*a - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-a + 1)*q^18 + (-3*a + 3)*q^19 + O(q^20)]]

    The number of coefficients stored can be changed::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(17); nfs
        [q - q^2 - q^4 - 2*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)[0].coefficient(79)
        Traceback (most recent call last):
        ...
        ValueError: The 79-th coefficient is not stored.
        sage: save_newforms(nfs, f, coefficients=100)
        sage: load_newforms(f)[0].coefficient(79)
        12
        sage: save_newforms(nfs, f, coefficients=[79])
        sage: load_newforms(f)[0].coefficient(79)
        12
        sage: load_newforms(f)[0].coefficient(78)
        Traceback (most recent call last):
        ...
        ValueError: The 78-th coefficient is not stored.

    Storing whether a curve has complex multiplication is optional::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(19); nfs
        [q - 2*q^3 - 2*q^4 + 3*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)[0].has_cm()
        False
        sage: save_newforms(nfs, f, save_cm=False)
        sage: load_newforms(f)[0].has_cm()
        Traceback (most recent call last):
        ...
        ValueError: Undetermined whether this newform has CM.

    Storing newforms also works for non-classical newforms::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: K = QuadraticField(3)
        sage: nfs = get_newforms(K.ideal(11), base_field=K, algorithm='magma')
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f, coefficients=15, save_cm=False)
        sage: nfs_copy = load_newforms(f)
        sage: nfs[0].base_field()
        Number Field in a with defining polynomial x^2 - 3
        sage: nfs_copy[0].base_field()
        Number Field in a with defining polynomial x^2 - 3
        sage: [nf.trace_of_frobenius(nf.base_field().prime_above(3)) for nf in nfs]
        [2, 1, -1, -2, -K1^2 + 5, K1^2 - 5]
        sage: [nf.trace_of_frobenius(nf.base_field().prime_above(3)) for nf in nfs_copy]
        [2, 1, -1, -2, -a^2 + 5, a^2 - 5]

    """
    if coefficients in ZZ and coefficients > Integer(0) :
        coefficients = list(range(coefficients))
    if repr_coefficients:
        coefficients = list(range(Integer(20) )) + [c for c in coefficients if c > Integer(20) ]
    else:
        coefficients = [c for c in coefficients]
    with open(file_name, "w+") as f:
        _write_element(newforms, f, coefficients, save_cm)

def _write_list(ls, f, coefficients, save_cm, indent=Integer(0) , indent_start=True):
    r"""Write a list to a file.

    Uses the rules specified in :func:`save_newforms` to write a list
    to a file.

    INPUT:

    - ``ls`` -- The list to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    - ``indent_start`` -- A boolean indicating whether the first
      symbol written should be indented.

    """
    if indent_start:
        f.write(" "*Integer(4) *indent)
    f.write('[')
    write_comma = False
    for element in ls:
        if write_comma:
            f.write(',')
        f.write('\n')
        _write_element(element, f, coefficients, save_cm, indent=indent+Integer(1) )
        write_comma=True
    f.write(']')

def _write_element(element, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write an element to a file.

    Uses the rules specified in :func:`save_newforms` to write an
    element to a file.

    INPUT:

    - ``element`` -- The element to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    if isinstance(element, WrappedNewform):
        _write_newform(element, f, coefficients, save_cm, indent=indent)
    elif is_DirichletCharacter(element):
        _write_character(element, f, coefficients, save_cm, indent=indent)
    elif is_Ideal(element):
        _write_ideal(element, f, coefficients, save_cm, indent=indent)
    elif is_NumberField(element) or element == QQ:
        _write_field(element, f, coefficients, save_cm, indent=indent)
    elif isinstance(element, LabeledElement):
        _write_labeled_element(element, f, coefficients, save_cm,
                               indent=indent)
    elif is_Polynomial(element):
        _write_polynomial(element, f, coefficients, save_cm, indent=indent)
    elif element in QQ or element in ZZ:
        _write_rational(element, f, indent=indent)
    elif hasattr(element, '__iter__'):
        _write_list(element, f, coefficients, save_cm, indent=indent)
    else:
        raise ValueError("Do not know how to write %s to file."%(element,))

def _write_newform(nf, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a newform to a file.

    Uses the rules specified in :func:`save_newforms` to write a
    newform to a file.

    INPUT:

    - ``nf`` -- The newform to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    if save_cm:
        cm = LabeledElement('cm', ZZ(nf.has_cm()))
    else:
        cm = LabeledElement('cm', ZZ(-Integer(1) ))
    character = nf.character()
    basefield = nf.base_field()
    if not basefield.is_absolute():
        basefield = basefield.absolute_field(names='a')
    level = LabeledElement('level', ZZ(nf.level()) if basefield == QQ else [nf.level()])
    coefficientfield = nf.coefficient_field()
    if not coefficientfield.is_absolute():
        coefficientfield = coefficientfield.absolute_field(names='a')
    if basefield == QQ:
        values = [(n, (QQ(nf.coefficient(n)) if nf.coefficient(n) in QQ
                       else LabeledElement('element',
                                           coefficientfield(nf.coefficient(n)).list())))
                  for n in coefficients]
    else:
        primes = [ZZ(p) for p in coefficients if ZZ(p).is_prime()]
        values = [(P, LabeledElement('element',
                                     coefficientfield(nf.trace_of_frobenius(P)).list()))
                  for p in primes for P in basefield.primes_above(p)
                  if not P.divides(nf.level())]
    basefield = LabeledElement('basefield', [basefield])
    coefficientfield = LabeledElement('coefficientfield', [coefficientfield])
    values = LabeledElement('values', values)
    element = LabeledElement('newform', [level, cm, character, basefield,
                                         coefficientfield, values])
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_character(eps, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a character to a file.

    Uses the rules specified in :func:`save_newforms` to write a
    character to a file.

    INPUT:

    - ``eps`` -- The character to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    eps = eps.primitive_character()
    conductor = LabeledElement('conductor', eps.conductor())
    ls_values = zip(eps.parent().unit_gens(), eps.element())
    values = LabeledElement('values', ls_values)
    element = LabeledElement('character', [conductor, values])
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_ideal(ideal, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a fractional ideal to a file.

    Uses the rules specified in :func:`save_newforms` to write a field
    to a file.

    INPUT:

    - ``ideal`` -- The fractional ideal to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    element = LabeledElement('ideal', [LabeledElement('element', g.list())
                                       for g in ideal.gens()])
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_field(field, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a field to a file.

    Uses the rules specified in :func:`save_newforms` to write a field
    to a file.

    INPUT:

    - ``field`` -- The field to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    if field is QQ:
        polynomial = PolynomialRing(QQ, names='x').gen()
    else:
        polynomial = field.defining_polynomial()
    element = LabeledElement('field', [polynomial])
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_polynomial(poly, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a polynomial to a file.

    Uses the rules specified in :func:`save_newforms` to write a
    polynomial to a file.

    INPUT:

    - ``poly`` -- The polynomial to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    - ``indent_start`` -- A boolean indicating whether the first
      symbol written should be indented.

    """
    element = LabeledElement('polynomial', poly.list())
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_rational(q, f, indent=Integer(0) , indent_start=True):
    r"""Write a rational to a file.

    Uses the rules specified in :func:`save_newforms` to write a
    rational to a file.

    INPUT:

    - ``q`` -- The rational to be written

    - ``f`` -- The file to be written to

    - ``indent`` -- An integer indicating the level of indentation to be used

    - ``indent_start`` -- A boolean indicating whether the first
      symbol written should be indented.

    """
    if indent_start:
        f.write(" "*Integer(4) *indent)
    f.write(str(QQ(q)))

def _write_labeled_element(element, f, coefficients, save_cm, indent=Integer(0) ):
    r"""Write a labeled element to a file.

    Uses the rules specified in :func:`save_newforms` to write a
    labeled element to a file.

    INPUT:

    - ``element`` -- The labeled element to be written

    - ``f`` -- The file to be written to

    - ``coefficients`` -- The indices of coefficients of newforms to
      be saved

    - ``save_cm`` -- A boolean indicating whether information about
      newforms having CM should be saved

    - ``indent`` -- An integer indicating the level of indentation to be used

    """
    f.write(" "*Integer(4) *indent)
    f.write("<%s> := "%(element.label,))
    if element.element in QQ:
        _write_rational(element.element, f, indent=indent, indent_start=False)
    elif hasattr(element.element, "__iter__"):
        _write_list(element.element, f, coefficients, save_cm, indent=indent,
                    indent_start=False)
    else:
        raise ValueError(str(elemnet.label) + ", " + str(element.element) +
                         " is not a valid labeled element")

def load_newforms(file_name):
    r"""Load newforms from a file.

    Load a newform or a list of newforms to a file. This file should
    contain for each newform information about its level, its
    character and some fourier coefficients. It can also contain
    whether or not the newform has complex multiplication, but this is
    optional.

    The file in which the newforms are stored should use the following
    notation, written here as regular expressions:

    <list> := '[' ( <element> ( ',' <element> )* )? ']'
    <element> := ( '<' <identifier> '>' ':=' )? ( <list> | <rational> )
    <identifier> := <letter>+
    <rational> := <integer> ( '/' <positive_integer> )?
    <integer> := ( '-' )? <zero> | <positive_integer>
    <positive_integer> := <non_zero_digit> ( <digit> )*
    <digit> := <zero> | <non_zero_digit>
    <non_zero_digit> := [1-9]
    <zero> := '0'
    <letter> := [a-zA-Z]

    Note that for <list>, <element> and <identifier> whitespace
    between the different building blocks is ignored. Furthermore it
    should represent different bits of data in the following way. Note
    that the function :func:`save_newforms` creates a file that
    satisfies all these properties.

    - A boolean value is represented by the integer 1 if it is True
      and the integer 0 if it is False.

    - An element of a number field is represented as the list of
      rational coefficients with respect to the power basis in the
      generator, preceded by the identifier 'element'

    - A polynomial with rational coefficients is represented by a list
      of its coefficients (starting at the constant term) preceded by
      the identifier 'polynomial'

    - A number field is represented by a list containing a polynomial,
      preceded by the identifier 'field'. The polynomial is the
      defining polynomial of the number field.

    - A (fractional) ideal of a number field is represented by a list
      containing the elements of the number field that generate the
      ideal, preceded by the identifier 'ideal'

    - A list of values of a function is represented by a list
      containing lists of exactly two elements, all preceded by the
      identifier 'values'. The function maps the first element of a
      list in the corresponding list to the second element thereof.

    - A character is represented by a list containing an integer
      preceded by the identifier 'conductor' and a list of values of
      the character preceded by the identifier 'values', all preceded
      by the identifier 'character'. The integer with identifier
      'conductor' will be the conductor of the character. The entry
      labeled values will be pairs of integers, such that if $\zeta$
      is the relevant $n$-th root of unity a pair $(k, e)$ appears in
      this list if the character takes the value $zeta^e$ at $k$.

    - A newform is represented by a list containing an integer or an
      ideal preceded by the identifier 'level', a boolean preceded by
      the identifier 'cm', a character, a number field preceded by the
      identifier 'basefield', a number field preceded by the
      identifier 'coefficientfield' and a list of values preceded by
      the identifier 'value', all preceded by the identifier
      'newform'. The first entry will be the level of the newform, the
      second a boolean indicating whether or not this newform has
      complex multiplication, the third the corresponding character,
      the fourth the base field of the newform, the fifth the
      coefficient field of the newform and the last the coefficients
      of the newform at (some) indices and the traces of frobenius at
      some prime ideals. For backwards compatibility the label
      basefield may be left out (in which case it is interpreted to be
      the rationals) and the label 'coefficientfield' may be replaced
      with an entry that is just a number field. The entry with label
      'cm' may be left out or set to -1 to indicate that this
      information is not known.

    - A list of things is represented as a list of the corresponding
      representations.

    INPUT:

    - ``file_name`` -- A string containing the file name from which
      the given newforms should be loaded.

    OUTPUT:

    An instance of :class:`WrappedNewform_stored` or a list thereof
    representing the newforms found in the given file.

    EXAMPLES::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(26); nfs
        [q - q^2 + q^3 + q^4 - 3*q^5 + O(q^6), q + q^2 - 3*q^3 + q^4 - q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)
        [q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + q^12 + q^13 + q^14 - 3*q^15 + q^16 - 3*q^17 + 2*q^18 + 2*q^19 + O(q^20),
         q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 - 3*q^12 - q^13 + q^14 + 3*q^15 + q^16 - 3*q^17 + 6*q^18 + 6*q^19 + O(q^20)]

    One can store multiple newform lists in a single file::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs1 = get_newforms(26); nfs1
        [q - q^2 + q^3 + q^4 - 3*q^5 + O(q^6), q + q^2 - 3*q^3 + q^4 - q^5 + O(q^6)]
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nfs2 = get_newforms(16, character=eps); nfs2
        [q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms([nfs1, nfs2], f)
        sage: load_newforms(f)
        [[q - q^2 + q^3 + q^4 - 3*q^5 - q^6 - q^7 - q^8 - 2*q^9 + 3*q^10 + 6*q^11 + q^12 + q^13 + q^14 - 3*q^15 + q^16 - 3*q^17 + 2*q^18 + 2*q^19 + O(q^20),
          q + q^2 - 3*q^3 + q^4 - q^5 - 3*q^6 + q^7 + q^8 + 6*q^9 - q^10 - 2*q^11 - 3*q^12 - q^13 + q^14 + 3*q^15 + q^16 - 3*q^17 + 6*q^18 + 6*q^19 + O(q^20)],
         [q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + (-2*a - 2)*q^12 + (a - 1)*q^13 + (2*a - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-a + 1)*q^18 + (-3*a + 3)*q^19 + O(q^20)]]

    The number of coefficients stored can be changed::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(17); nfs
        [q - q^2 - q^4 - 2*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)[0].coefficient(79)
        Traceback (most recent call last):
        ...
        ValueError: The 79-th coefficient is not stored.
        sage: save_newforms(nfs, f, coefficients=100)
        sage: load_newforms(f)[0].coefficient(79)
        12
        sage: save_newforms(nfs, f, coefficients=[79])
        sage: load_newforms(f)[0].coefficient(79)
        12
        sage: load_newforms(f)[0].coefficient(78)
        Traceback (most recent call last):
        ...
        ValueError: The 78-th coefficient is not stored.

    Storing whether a curve has complex multiplication is optional::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: nfs = get_newforms(19); nfs
        [q - 2*q^3 - 2*q^4 + 3*q^5 + O(q^6)]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(nfs, f)
        sage: load_newforms(f)[0].has_cm()
        False
        sage: save_newforms(nfs, f, save_cm=False)
        sage: load_newforms(f)[0].has_cm()
        Traceback (most recent call last):
        ...
        ValueError: Undetermined whether this newform has CM.

    """
    with open(file_name, 'r') as f:
        return _interpret_element(_read_element(f))

def _interpret_element(element, field=None):
    r"""Recover a data structure from its file representation.

    The way data is represented can be found in the description of
    :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    - `field` -- A number field relevant to interpreting this element.

    OUTPUT:

    The data structure that the given element represents.

    """
    if isinstance(element, LabeledElement):
        return _interpret_labeled_element(element, field=field)
    elif isinstance(element, list):
        return [_interpret_element(subelement, field=field)
                for subelement in element]
    else:
        return element

def _interpret_labeled_element(element, field=None):
    r"""Recover a data structure from its file representation.

    The way data is represented can be found in the description of
    :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    - `field` -- A number field relevant to interpreting this element.

    OUTPUT:

    The data structure that the given labeled element represents.

    """
    if element.label == 'newform':
        return _interpret_newform(element.element)
    elif element.label == 'character':
        return _interpret_character(element.element)
    elif element.label == 'element':
        return _interpret_field_element(element.element, field=field)
    elif element.label == 'ideal':
        return _interpret_ideal(element.element, field=field)
    elif element.label == 'field':
        return _interpret_field(element.element)
    elif element.label == 'polynomial':
        return _interpret_polynomial(element.element)
    else:
        return element

def _interpret_polynomial(element):
    r"""Recover a polynomial from its file representation.

    The way polynomials are represented can be found in the
    description of :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    OUTPUT:

    The polynomial over $\QQ$ that the given labeled element
    represents.

    """
    R = QQ['x']; (x,) = R._first_ngens(1)
    return R(element)

def _interpret_field(element):
    r"""Recover a field from its file representation.

    The way fields are represented can be found in the description of
    :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    OUTPUT:

    The number field that the given labeled element represents.

    """
    if not isinstance(element, list):
        raise ValueError("%s is not a list."%(element,))
    if len(element) != Integer(1) :
        raise ValueError("%s does not have length 1."%(element,))
    polynomial = _interpret_element(element[Integer(0) ])
    if not is_Polynomial(polynomial):
        raise ValueError("%s is not a polynomial."%(polynomial,))
    if polynomial.degree() == Integer(1) :
        return QQ
    return NumberField(polynomial, names='a')

def _interpret_ideal(element, field=None):
    r"""Recover a fractional ideal from its file representation.

    The way fractional ideals are represented can be found in the
    description of :func:`load_newforms`.

    INPUT:

    - `element` -- An element

    - `field` -- A number field in which this ideal would reside

    OUTPUT:

    The fractional ideal that the given labeled element represents.

    """
    if not isinstance(element, list):
        raise ValueError("%s is not a list."%(element,))
    if field is None:
        raise ValueError("Should specify a field to " +
                         "interpret the ideal %s."%(element,))
    gens = [_interpret_element(gen, field=field) for gen in element]
    return field.ideal(gens)

def _interpret_field_element(element, field=None):
    r"""Recover a field element from its file representation.

    The way field elements are represented can be found in the
    description of :func:`load_newforms`.

    INPUT:

    - `element` -- An element

    - `field` -- A number field in which this ideal would reside

    OUTPUT:

    The field element that the given labeled element represents.

    """
    if not isinstance(element, list) or not all(a in QQ for a in element):
        raise ValueError("%s is not a list of rationals."%(element,))
    if field is None:
        raise ValueError("Should specify a field to " +
                         "interpret the field element %s."%(element,))
    return field(element)

def _interpret_character(element):
    r"""Recover a character from its file representation.

    The way characters are represented can be found in the
    description of :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    OUTPUT:

    The Dirichlet character that the given labeled element represents.

    """
    if not isinstance(element, list):
        raise valueError("%s is not a list."%(element,))
    conductor=None
    values=None
    for part in element:
        if isinstance(part, LabeledElement):
            if (part.label.lower() == 'conductor' and conductor is None and
                part.element in ZZ):
                conductor = ZZ(part.element)
            elif (part.label.lower() == 'values' and values is None and
                  isinstance(part.element, list)):
                values=dict()
                for pair in part.element:
                    if (not isinstance(pair, list) or len(pair) != Integer(2)  or
                        pair[Integer(0) ] not in ZZ):
                        raise ValueError("Expected a pair for character " +
                                         "values, but got " + str(pair))
                    if pair[Integer(1) ] in ZZ:
                        values[ZZ(pair[Integer(0) ])] = ZZ(pair[Integer(1) ])
                    else:
                        raise ValueError("Expected a pair for character " +
                                         "values, but got " + str(pair))
            else:
                raise ValueError("Unexpected element " + str(part.element) +
                                 " with label " + str(part.label) +
                                 " for character.")
        else:
            raise ValueError("Unexpected element %s for character."%(part,))
    if conductor is None or values is None:
        raise ValueError("Not enough arguments to make a character.")
    D = DirichletGroup(conductor)
    try:
        pows = D._zeta_powers
        return D([pows[values[n]] for n in D.unit_gens()])
    except KeyError as e:
        raise ValueError("Requires value at " + str(e) +
                         " to construct Dirichlet character.")

def _interpret_newform(element):
    r"""Recover a newform from its file representation.

    The way newforms are represented can be found in the
    description of :func:`load_newforms`.

    INPUT:

    - `element` -- A labeled element

    OUTPUT:

    The newform that the given labeled element represents.

    """
    if not isinstance(element, list):
        raise valueError("%s is not a list."%(element,))
    level=None
    cm=None
    character=None
    basefield=None
    coefficientfield=None
    values=None
    for part in element:
        if isinstance(part, LabeledElement):
            if (part.label.lower() == 'character' and character is None):
                character = _interpret_character(part.element)
            elif (part.label.lower() == 'field' and coefficientfield is None):
                coefficientfield = _interpret_field(part.element)
            elif (part.label.lower() == 'coefficientfield' and coefficientfield is None
                  and isinstance(part.element, list) and len(part.element) == Integer(1) ):
                coefficientfield = _interpret_element(part.element[Integer(0) ])
            elif (part.label.lower() == 'basefield' and basefield is None
                  and isinstance(part.element, list) and len(part.element) == Integer(1) ):
                basefield = _interpret_element(part.element[Integer(0) ])
            elif (part.label.lower() == 'level' and level is None):
                level = part.element
            elif (part.label.lower() == 'values' and values is None and
                  isinstance(part.element, list)):
                values = part.element
            elif (part.label.lower() == 'cm' and cm is None and
                  part.element in ZZ):
                if part.element != -Integer(1) : # -1 means cm is undefined
                    cm = (part.element != Integer(0) )
            else:
                raise ValueError("Unexpected element " + str(part.element) +
                                 " with label " + str(part.label) +
                                 " for newform.")
        else:
            raise ValueError("Unexpected element %s for newform."%(part,))
    if basefield is None:
        basefield = QQ # Included for backwards compatibility
    if (level is None or character is None or basefield is None or coefficientfield is None):
        raise ValueError("Not enough arguments to make a newform.")
    if basefield != QQ:
        level = _interpret_element(level[Integer(0) ], field=basefield)
    if not (level in ZZ or (is_Ideal(level) and level in basefield.ideal_monoid())):
        raise ValueError("%s can not be the level of a newform"%(level,))
    coefficients = {}
    traces = {}
    for pair in values:
        if (not isinstance(pair, list) or len(pair) != Integer(2) ):
            raise ValueError("Expected a pair for newform " +
                             "coefficients, but got " + str(pair))
        arg = _interpret_element(pair[Integer(0) ], field=basefield)
        val = _interpret_element(pair[Integer(1) ], field=coefficientfield)
        if not val in coefficientfield:
            raise ValueError("Expectend an element of %s"%(coefficientfield,)
                             + ", but got %s instead"%(val,))
        if arg in ZZ:
            coefficients[ZZ(arg)] = val
        elif is_Ideal(arg) and arg in basefield.ideal_monoid():
            traces[arg] = val
        else:
            raise ValueError("Expected an integer or an element of " +
                             str(basefield) + ", but got %s instead"%(arg))
    return WrappedNewform_stored(basefield, level, character,
                                 coefficientfield, coefficients, traces, cm)

def _read_element(f):
    r"""Read an element from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The element read from from the file.

    """
    whitespace = re.compile('\\s')
    label = None
    element = None
    s = f.read(Integer(1) )
    while len(s) > Integer(0) :
        if whitespace.match(s):
            pass
        elif s == '<':
            f.seek(f.tell()-Integer(1) )
            label = _read_identifier(f)
            _read_colon_equals(f)
        elif s == '[':
            f.seek(f.tell()-Integer(1) )
            element = _read_list(f)
            break
        elif s.isdigit() or s == '-':
            f.seek(f.tell()-Integer(1) )
            element = _read_rational(f)
            break
        else:
            raise ValueError("Read unexpected character " + str(s) +
                             " while reading an element.")
        s = f.read(Integer(1) )
    else:
        raise ValueError("File ended before reading an element.")
    if label is None:
        return element
    else:
        return LabeledElement(label, element)

def _read_colon_equals(f):
    r"""Read ':=' from a file.

    INPUT:

    - ``f`` -- The file to be read from

    """
    whitespace = re.compile('\\s')
    s = f.read(Integer(1) )
    while len(s) > Integer(0) :
        if whitespace.match(s):
            pass
        elif s == ':':
           s = f.read(Integer(1) )
           if s == '=':
               return
           else:
               raise ValueError("Attempting to read ':=' at '" + str(s) +
                                "', but failed.")
        else:
           raise ValueError("Attempting to read ':=' at '" + str(s) +
                            "', but failed.")
        s = f.read(Integer(1) )
    raise ValueError("Attempting to read ':=', but encountered end of file.")

def _read_identifier(f):
    r"""Read an identifier from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The string that is the identifier read.

    """
    s = f.read(Integer(1) )
    if s != '<':
        raise ValueError("Attempting to read '<', but read %s"%(s,))
    result = ""
    s = f.read(Integer(1) )
    while len(s) > Integer(0) :
        if s.isalpha():
            result = result + s
        elif s == '>':
            return result
        else:
            raise ValueError("Attempting to read a letter, but read %s"%(s,))
        s = f.read(Integer(1) )
    raise ValueError("Reached end of file whilst reading an identifier.")

def _read_list(f):
    r"""Read a list from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The list read from the file.

    """
    whitespace = re.compile('\\s')
    s = f.read(Integer(1) )
    if s != '[':
        raise ValueError("Attempting to read the start of a list at '" +
                         str(s) + "', but failed")
    result = []
    s = f.read(Integer(1) )
    can_close = True
    need_comma = False
    while len(s) > Integer(0) :
        if whitespace.match(s):
            pass
        elif s == ']' and can_close:
            return result
        elif s == ',' and need_comma:
            need_comma = False
            can_close = False
        elif not need_comma:
            f.seek(f.tell()-Integer(1) )
            result.append(_read_element(f))
            need_comma = True
            can_close = True
        else:
            raise ValueError("Attempting to read more of a list, " +
                             "but encountered '" + str(s) + "'")
        s = f.read(Integer(1) )
    raise ValueError("Encountered end of file whilst reading a list.")

def _read_rational(f):
    r"""Read a rational number from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The rational number read from the file.

    """
    result = _read_integer(f)
    s = f.read(Integer(1) )
    if s == '/':
        result = result + '/' + _read_positive_integer(f)
    else:
        f.seek(f.tell()-Integer(1) )
    return QQ(result)

def _read_integer(f):
    r"""Read an integer from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The integer read from the file as a string.

    """
    result = ''
    s = f.read(Integer(1) )
    if s == '-':
        result = result + '-'
    else:
        f.seek(f.tell()-Integer(1) )
    try:
        result = result + _read_zero(f)
    except ValueError:
        f.seek(f.tell()-Integer(1) )
        result = result + _read_positive_integer(f)
    return result

def _read_positive_integer(f):
    r"""Read a positive integer from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The positive integer read from the file as a string.

    """
    result = _read_non_zero_digit(f)
    try:
        while True:
            result = result + _read_digit(f)
    except ValueError:
        f.seek(f.tell()-Integer(1) )
        return result

def _read_digit(f):
    r"""Read a digit from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The digit read from the file as a string.

    """
    s = f.read(Integer(1) )
    if not s.isdigit():
        raise ValueError("Attempting to read a digit, but read %s"%(s,))
    return s

def _read_non_zero_digit(f):
    r"""Read a non-zero digit from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The non-zero digit read from the file as a string.

    """
    s = f.read(Integer(1) )
    if s.isdigit() and s != Integer(0) :
        return s
    else:
        raise ValueError("Attempting to read a non-zero digit, but read %s"%(s,))

def _read_zero(f):
    r"""Read the digit '0' from a file.

    The way data is stored in and hence read from the file can be
    found in the description of :func:`load_newforms`.

    INPUT:

    - ``f`` -- The file to be read from

    OUTPUT:

    The string '0'

    """
    s = f.read(Integer(1) )
    if s != '0':
        raise ValueError("Attempting to read 0, but read %s"%(s,))
    return s

class WrappedNewform(SageObject):
    r"""A wrapper class around a newform of weight 2.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    This class is a base class, but should not be used as an
    instance. It rather provides a template for all classes that
    inherit it.

    EXAMPLES:

    A wrapper around a Sage newform::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nf = get_newforms(16, character=eps)[0]; nf
        q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    A wrapper around a magma newform::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nf = get_newforms(16, character=eps, algorithm='magma')[0]; nf
        q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + O(q^12)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    A wrapper around a newform from a file::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(get_newforms(16, character=eps), f)
        sage: nf = load_newforms(f)[0]; nf
        q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + (-2*a - 2)*q^12 + (a - 1)*q^13 + (2*a - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-a + 1)*q^18 + (-3*a + 3)*q^19 + O(q^20)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    """

    def __init__(self):
        r""" Initialize this object"""
        self._embeddings = {}

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        A non-negative integer describing the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        raise NotImplementedError()

    def can_compute_frobenius(self, prime):
        r"""Determine whether the trace of Frobenius can be computed for `prime`.

        INPUT:

        - ``prime`` -- A finite prime of the :meth:`base_field`. If
          :meth:`base_field` is the rationals it should be a prime
          number, otherwise a prime ideal.

        OUTPUT:

        - `True`, if this wrapped newform can compute the trace of
          frobenius at the prime `prime`. At least the `prime` should
          then not divide the :meth:`level`.

        - `False`, otherwise.

        """
        return not prime.divides(self.level())

    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        raise NotImplementedError()

    def base_field(self):
        r"""Give the base field for this newform.

        For classical newforms this is always the rationals. For
        Hilbert modular forms this is the totally real field for which
        this is a Hilbert modular form. For Bianchi modular forms it
        is the complex field for which this is a Bianchi modular form.

        OUTPUT:

        The rational field if this newform is a classical modular form.

        The totally real field for which this newform is a Hilbert
        modular form if it is a Hilbert modular form.

        The complex field for which this newform is a Bianchi modular
        form if it is a Bianchi modular form.

        """
        return QQ

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        raise NotImplementedError()

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the eigenvalues of the Hecke operators on
        this newform are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient_field()
            Rational Field
            sage: nf = get_newforms(31)[0]
            sage: nf.coefficient_field()
            Number Field in a0 with defining polynomial x^2 - x - 1

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        raise NotImplementedError()

    def set_embedding(self, field, embedding=None):
        r"""Set the embedding of the coefficient field into a given field

        Sets the embedding returned by :meth:`embedding` for a given
        field to a given embedding. If no embedding is given it will
        just set it to a embedding. If no embedding exists, will raise
        a ValueError.

        INPUT:

        - ``field`` -- A number field

        - ``embedding`` -- An embedding of the coefficient field of
          this newform into the given field `field` or None
          (default). If set to None, will use the method embedding to
          find an embedding of the coefficient field into the given
          field or raise a ValueError if none exist.

        EXAMPLES::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: f = get_newforms(2^3*3^4)[4]
            sage: K = f.coefficient_field(); K
            Number Field in a4 with defining polynomial x^2 + 8*x + 4
            sage: L = CyclotomicField(12)
            sage: f.embedding(L)
            Ring morphism:
              From: Number Field in a4 with defining polynomial x^2 + 8*x + 4
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: a4 |--> 2*zeta12^3 - 4*zeta12 - 4
            sage: f.set_embedding(L, K.embeddings(L)[1])
            sage: f.embedding(L)
            Ring morphism:
              From: Number Field in a4 with defining polynomial x^2 + 8*x + 4
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: a4 |--> -2*zeta12^3 + 4*zeta12 - 4
            sage: f.set_embedding(QuadraticField(5))
            Traceback (most recent call last):
            ...
            ValueError: No embedding from Number Field in a4 with defining polynomial x^2 + 8*x + 4 to Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? exists

        """
        if embedding is None:
            embeddings = self.coefficient_field().embeddings(field)
            if len(embeddings) == Integer(0) :
                raise ValueError("No embedding from " +
                                 str(self.coefficient_field()) +
                                 " to " + str(field) + " exists")
            embedding = embeddings[Integer(0) ]
        self._embeddings[field] = embedding

    def embedding(self, field):
        r"""Give an embedding of the coefficient field into the given field

        The given embedding is always the same and can be set using
        the method :meth:`set_embedding`.

        INPUT:

        - ``field`` -- A number field.

        OUTPUT:

        An embedding of the coefficient field of this newform into the
        given field. Raises a ValueError if no such embedding exists.

        .. SEEALSO:

            :meth:`set_embedding`
            :meth:`coefficient_field`

        EXAMPLES::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: f = get_newforms(2^3*3^4)[4]
            sage: K = f.coefficient_field(); K
            Number Field in a4 with defining polynomial x^2 + 8*x + 4
            sage: L = CyclotomicField(12)
            sage: f.embedding(L)
            Ring morphism:
              From: Number Field in a4 with defining polynomial x^2 + 8*x + 4
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: a4 |--> 2*zeta12^3 - 4*zeta12 - 4
            sage: f.set_embedding(L, K.embeddings(L)[1])
            sage: f.embedding(L)
            Ring morphism:
              From: Number Field in a4 with defining polynomial x^2 + 8*x + 4
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: a4 |--> -2*zeta12^3 + 4*zeta12 - 4
            sage: f.embedding(QuadraticField(5))
            Traceback (most recent call last):
            ...
            ValueError: No embedding from Number Field in a4 with defining polynomial x^2 + 8*x + 4 to Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? exists

        """
        if not (field in self._embeddings):
            self.set_embedding(field)
        return self._embeddings[field]

    def q_expansion(self, prec=Integer(20) ):
        """Give the q-expansion of this newform.

        INPUT:

        - ``prec`` -- A non-negative integer (default: 20) giving a
          bound on the powers that should be present in this
          q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given as a power
        series in q with coefficients in the coefficient field of this
        newform and capped at the given precision `prec`.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.q_expansion()
            q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + 3*q^11 + 4*q^12 - 4*q^13 - 6*q^15 + 4*q^16 - 3*q^17 + q^19 + O(q^20)
            sage: nf.q_expansion(10)
            q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + O(q^10)

        .. SEE_ALSO::

            :meth:`coefficient`
            :meth:`coefficient_field`

        """
        R = self.coefficient_field()[['q']]; (q,) = R._first_ngens(1)
        result = sum(self.coefficient(n) * q**n for n in range(prec))
        return result.add_bigoh(prec)

    @cached_method
    def _trace_power_formula(self, power):
        r"""Give the formula to compute the trace of frobenius to a given
        power.

        The trace of the galois representation of this newform at a
        frobenius element to the power $n$ can be computed from the
        trace and determinant at the frobenius element with this
        formula

        INPUT:

        - ``power`` -- The power $n$ of the frobenius element.

        """
        R = QQ['x, y']; (x, y,) = R._first_ngens(2)
        return polynomial_to_symmetric(x**power + y**power)

    def trace_of_frobenius(self, prime, power=Integer(1) ):
        """Give the trace of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The trace of $\rho(F_p^n)$, where $\rho$ is the mod $l$ or
        l-adic galois representation associated to this newform, $F_p$
        is the frobenius element at the given prime, and $n$ is the
        given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.trace_of_frobenius(2)
            0
            sage: nf.trace_of_frobenius(7)
            -1
            sage: nf.trace_of_frobenius(7, power=2)
            -13

        .. SEE_ALSO::

            :meth:`determinant_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == Integer(1) :
            return self.coefficient(prime)
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=Integer(1) ):
        """Give the determinant of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The determinant of $\rho(F_p^n)$, where $\rho$ is the mod $l$
        or $l$-adic galois representation associated to this newform,
        $F_p$ is the frobenius element at the given prime, and $n$ is
        the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.determinant_of_frobenius(5)
            5
            sage: nf.determinant_of_frobenius(7)
            7
            sage: nf.determinant_of_frobenius(7, power=2)
            49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        D = self.character()(prime) * prime
        return D**power

    def characteristic_polynomial(self, prime, power=Integer(1) ):
        """Give the characteristic polynomial of the frobenius element acting
        on this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The characteristic polynomial of $\rho(F_p^n)$, where $\rho$
        is the mod $l$ or $l$-adic galois representation associated to
        this newform, $F_p$ is the frobenius element at the given
        prime, and $n$ is the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.characteristic_polynomial(5)
            x^2 - 3*x + 5
            sage: nf.characteristic_polynomial(7)
            x^2 + x + 7
            sage: nf.characteristic_polynomial(7, power=2)
            x^2 + 13*x + 49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`determinant_of_frobenius`

        """
        T = self.trace_of_frobenius(prime, power=power)
        D = self.determinant_of_frobenius(prime, power=power)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        R = K['x']; (x,) = R._first_ngens(1)
        return x**Integer(2)  - T_map(T)*x + D_map(D)

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        raise NotImplementedError()

    qexp = q_expansion

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform()

    def twist(self, character):
        r"""Twist this newform with a character.

        INPUT:

        - ``character`` -- A Dirichlet character.

        OUTPUT:

        A newform that is this newform twisted by `character`.

        """
        if character.conductor() == Integer(1) :
            return self
        else:
            return WrappedNewform_twisted(newform=self, twist=character)

    def _repr_(self):

        """Give a string representation of this newform"""
        return str(self.q_expansion())

    def _latex_(self):

        return latex(self.q_expansion())

class WrappedNewform_sage(WrappedNewform):
    r"""A wrapper class around a Sage newform of weight 2.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    EXAMPLE::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nf = get_newforms(16, character=eps)[0]; nf
        q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + O(q^6)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    """

    def __init__(self, newform):
        r"""Initialize a wrapped newform.

        INPUT:

        - ``newform`` -- The sage newform that should be wrapped.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import WrappedNewform_sage
            sage: WrappedNewform_sage(Newforms(19)[0])
            q - 2*q^3 - 2*q^4 + 3*q^5 + O(q^6)

        """
        self._f = newform
        WrappedNewform.__init__(self)

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        A non-negative integer describing the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        return self._f.level()

    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        return self._f.character()

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        if n == Integer(0) :
            return self.coefficient_field()(Integer(0) )
        else:
            return self._f.coefficient(n)

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the coefficients of the q-expansion of
        this newform at infinity are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient_field()
            Rational Field
            sage: nf = get_newforms(31)[0]
            sage: nf.coefficient_field()
            Number Field in a0 with defining polynomial x^2 - x - 1

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._f.base_ring()

    def q_expansion(self, prec=Integer(20) ):
        """Give the q-expansion of this newform.

        INPUT:

        - ``prec`` -- A non-negative integer (default: 20) giving a
          bound on the powers that should be present in this
          q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given as a power
        series in q with coefficients in the coefficient field of this
        newform and capped at the given precision `prec`.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.q_expansion()
            q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + 3*q^11 + 4*q^12 - 4*q^13 - 6*q^15 + 4*q^16 - 3*q^17 + q^19 + O(q^20)
            sage: nf.q_expansion(10)
            q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + O(q^10)

        .. SEE_ALSO::

            :meth:`coefficient`
            :meth:`coefficient_field`

        """
        return self._f.q_expansion(prec=prec)

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        return self._f.has_cm()

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_sage(self._f)

    def _repr_(self):
        """Give a string representation of this newform."""
        return str(self._f)

    def _latex_(self):
        """Give a latex representation of this newform."""
        return latex(self._f)

    qexp = q_expansion

class WrappedNewform_magma(WrappedNewform):
    r"""A wrapper class around a magma newform of weight 2.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    EXAMPLE::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: nf = get_newforms(16, character=eps, algorithm='magma')[0]; nf
        q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + O(q^12)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    """

    def __init__(self, newform):
        r"""Initialize a wrapped newform

        INPUT:

        - ``newform`` -- A newform as a magma object

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import WrappedNewform_magma
            sage: cfs = magma.CuspForms(19)
            sage: WrappedNewform_magma(magma.Newforms(cfs)[1][1])
            q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + 3*q^11 + O(q^12)

        """
        self._f = newform
        WrappedNewform.__init__(self)

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        A non-negative integer describing the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        return self._f.Level().sage()

    @cached_method
    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        eps_f = self._f.DirichletCharacter()
        N = eps_f.Modulus().sage()
        N0 = eps_f.Conductor().sage()
        L = eps_f.BaseRing().sage()
        gens = Integers(N).unit_gens()
        for eps in DirichletGroup(N0, base_ring=L):
            if all(eps(n) == eps_f(n).sage() for n in gens):
                return eps
        raise ValueError("No sage character corresponds to %s"%(eps_f,))

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        return self._f.Coefficient(n).sage()

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the coefficients of the q-expansion of
        this newform at infinity are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient_field()
            Rational Field
            sage: nf = get_newforms(31)[0]
            sage: nf.coefficient_field()
            Number Field in a0 with defining polynomial x^2 - x - 1

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._f.BaseField().sage()

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        mds = self._f.ModularSymbols().NewformDecomposition()[Integer(1) ]
        return mds.HasCM(Proof=proof)

    def _repr_(self):
        """Give a string representation of this newform"""
        return str(self._f)

    def _latex_(self):
        """Give a latex representation of this newform."""
        return latex(self._f)

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_magma(self._f)

class WrappedNewform_stored(WrappedNewform):
    r"""A wrapper class around a newform of (parallel) weight 2 defined by
    stored data.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    The data that has to be provided in order to construct such a
    newform is the base field, the level, the character, the
    coefficient field. For functionality one should also provide some
    coefficients of the q-expansion for classical modular forms or
    some values for traces of Frobenius at specific primes. Optionally
    whether or not the newform has complex multiplication can also be
    provided.

    EXAMPLE::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms, save_newforms, load_newforms
        sage: eps = DirichletGroup(16).gens()[1]
        sage: f = tmp_filename(ext='.nfs')
        sage: save_newforms(get_newforms(16, character=eps), f)
        sage: nf = load_newforms(f)[0]; nf
        q + (-a - 1)*q^2 + (a - 1)*q^3 + 2*a*q^4 + (-a - 1)*q^5 + 2*q^6 - 2*a*q^7 + (-2*a + 2)*q^8 + a*q^9 + 2*a*q^10 + (a + 1)*q^11 + (-2*a - 2)*q^12 + (a - 1)*q^13 + (2*a - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-a + 1)*q^18 + (-3*a + 3)*q^19 + O(q^20)
        sage: nf.level()
        16
        sage: nf.character()
        Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

    """

    def __init__(self, base_field, level, character,
                 coefficient_field, coefficients={}, traces={},
                 cm=None):
        r"""Initialize a wrapped newform

        INPUT:

        - ``base_field`` -- The number field over which this is a
          newform. This would be the rationals for a classical modular
          form, a totally real field for a Hilbert modular
          form, or a complex field for a Bianchi modular form.

        - ``level`` -- A non-negative integer which is the level of
          the newform.

        - ``character`` -- A Dirichlet character which is the
          character of the newform.

        - ``coefficient_field`` -- The number field over which the
          q-expansion of this newform is defined.

        - ``coefficients`` -- A dictionary indexed by non-negative
          integers and with as values the corresponding fourier
          coefficients of the q-expansion of this newform at
          infinity. By default this is the empty dictionary.

        - ``traces`` -- A dictionary indexed by prime ideals of the
          base field and with as values the corresponding traces of
          Frobenius at those primes. By default this is the empty
          dictionary.

        - ``cm`` -- A boolean value or 'None' (default) indicating
          whether this newform has complex multiplication or whether
          this is unknown ('None').

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms, WrappedNewform_stored
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: K = nf.coefficient_field()
            sage: c = {n : nf.coefficient(n) for n in range(50)}
            sage: WrappedNewform_stored(QQ, 16, eps, K, coefficients=c, cm=nf.has_cm())
            q + (-zeta4 - 1)*q^2 + (zeta4 - 1)*q^3 + 2*zeta4*q^4 + (-zeta4 - 1)*q^5 + 2*q^6 - 2*zeta4*q^7 + (-2*zeta4 + 2)*q^8 + zeta4*q^9 + 2*zeta4*q^10 + (zeta4 + 1)*q^11 + (-2*zeta4 - 2)*q^12 + (zeta4 - 1)*q^13 + (2*zeta4 - 2)*q^14 + 2*q^15 - 4*q^16 - 2*q^17 + (-zeta4 + 1)*q^18 + (-3*zeta4 + 3)*q^19 + O(q^20)

        """
        self._base = base_field
        self._level = level
        self._eps = character
        self._K = coefficient_field
        self._coeffs = coefficients
        self._traces = traces
        self._cm = cm
        WrappedNewform.__init__(self)

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        A non-negative integer describing the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        return self._level

    def can_compute_frobenius(self, prime):
        r"""Determine whether the trace of Frobenius can be computed for `prime`.

        INPUT:

        - ``prime`` -- A finite prime of the :meth:`base_field`. If
          :meth:`base_field` is the rationals it should be a prime
          number, otherwise a prime ideal.

        OUTPUT:

        - `True`, if this wrapped newform can compute the trace of
          frobenius at the prime `prime`. At least the `prime` should
          then not divide the :meth:`level`.

        - `False`, otherwise.

        """
        return prime in self._coeffs or prime in self._traces

    @cached_method
    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        return self._eps

    def base_field(self):
        r"""Give the base field for this newform.

        For classical newforms this is always the rationals. For
        Hilbert modular forms this is the totally real field for which
        this is a Hilbert modular form. For Bianchi modular forms it
        is the complex field for which this is a Bianchi modular form.

        OUTPUT:

        The rational field if this newform is a classical modular form.

        The totally real field for which this newform is a Hilbert
        modular form if it is a Hilbert modular form.

        The complex field for which this newform is a Bianchi modular
        form if it is a Bianchi modular form.

        """
        return self._base

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        try:
            return self._coeffs[n]
        except KeyError:
            raise ValueError("The %s-th coefficient is not stored."%(n,))

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the coefficients of the q-expansion of
        this newform at infinity are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient_field()
            Rational Field
            sage: nf = get_newforms(31)[0]
            sage: nf.coefficient_field()
            Number Field in a0 with defining polynomial x^2 - x - 1

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._K

    def trace_of_frobenius(self, prime, power=Integer(1) ):
        """Give the trace of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The trace of $\rho(F_p^n)$, where $\rho$ is the mod $l$ or
        l-adic galois representation associated to this newform, $F_p$
        is the frobenius element at the given prime, and $n$ is the
        given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.trace_of_frobenius(2)
            0
            sage: nf.trace_of_frobenius(7)
            -1
            sage: nf.trace_of_frobenius(7, power=2)
            -13

        .. SEE_ALSO::

            :meth:`determinant_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        F = self.base_field()
        if F != QQ and not is_Ideal(prime):
            prime = F.ideal(prime)
        if not(prime in ZZ or prime in F.ideal_monoid()) or not prime.is_prime():
            raise ValueError("%s is not a valid prime."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == Integer(1) :
            if F == QQ:
                return self.coefficient(prime)
            else:
                try:
                    return self._traces[prime]
                except KeyError:
                    raise ValueError("The trace of Frobenius at %s is not stored."%(prime,))
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=Integer(1) ):
        """Give the determinant of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The determinant of $\rho(F_p^n)$, where $\rho$ is the mod $l$
        or $l$-adic galois representation associated to this newform,
        $F_p$ is the frobenius element at the given prime, and $n$ is
        the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.determinant_of_frobenius(5)
            5
            sage: nf.determinant_of_frobenius(7)
            7
            sage: nf.determinant_of_frobenius(7, power=2)
            49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        if self.base_field() == QQ:
            if prime not in ZZ or not prime.is_prime():
                raise ValueError("%s is not a prime number."%(prime,))
            if prime.divides(self.level()):
                raise ValueError("%s divides the level: %s."%(prime, self.level()))
            D = self.character()(prime) * prime
            return D**power
        raise NotImplementedError()

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        if self._cm is None:
            raise ValueError("Undetermined whether this newform has CM.")
        return self._cm

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_stored(self._base, self._level,
                                     self._eps, self._K, self._coeffs,
                                     self._traces, self._cm)

    def _repr_(self):
        """Give a string representation of this newform"""
        try:
            return str(self.q_expansion())
        except ValueError:
            return "Loaded newform with limited coefficients."

    def _latex_(self):
        """Give a latex representation of this newform."""
        try:
            return latex(self.q_expansion())
        except ValueError:
            return "\\text{Loaded newform with limited coefficients}"

class WrappedNewform_magma_hilbert(WrappedNewform):
    r"""A wrapper class around a magma Hilbert modular newform of parallel
    weight 2.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    EXAMPLE::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: K = QuadraticField(2)
        sage: nf = get_newforms(K.ideal(11), base_field=K, algorithm='magma')[0]; nf
        Element of Cuspidal newform space of Hilbert modular forms
        sage: nf.level()
        Fractional ideal (11)
        sage: nf.base_field()
        Number Field in a with defining polynomial x^2 - 2

    """

    def __init__(self, space):
        r"""Initialize a wrapped newform

        INPUT:

        - ``newform`` -- A Hilbert modular newform as a magma object

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import WrappedNewform_magma_hilbert
            sage: K = QuadraticField(2)
            sage: cfs = magma.HilbertCuspForms(K, K.ideal(9))
            sage: nfs = cfs.NewSubspace()
            sage: WrappedNewform_magma_hilbert(nfs.NewformDecomposition()[1])
            Element of Cuspidal newform space of Hilbert modular forms

        """
        self._f = space.Eigenform()
        self._M = space
        WrappedNewform.__init__(self)

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        An ideal of the totally real number field associated to this
        Hilbert modular form that describes the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        return self._M.Level().sage()

    @cached_method
    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        return DirichletGroup(Integer(1) )[Integer(0) ]

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        raise NotImplementedError();

    def base_field(self):
        r"""Give the base field for this newform.

        For classical newforms this is always the rationals. For
        Hilbert modular forms this is the totally real field for which
        this is a Hilbert modular form. For Bianchi modular forms it
        is the complex field for which this is a Bianchi modular form.

        OUTPUT:

        The rational field if this newform is a classical modular form.

        The totally real field for which this newform is a Hilbert
        modular form if it is a Hilbert modular form.

        The complex field for which this newform is a Bianchi modular
        form if it is a Bianchi modular form.

        """
        return self._M.BaseField().sage()

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the coefficients of the q-expansion of
        this newform at infinity are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: K = QuadraticField(2)
            sage: nf = get_newforms(K.ideal(5), base_field=K, algorithm='magma')[0]
            sage: nf.coefficient_field()
            Number Field in K1 with defining polynomial x^2 + 2*x - 2

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._M.HeckeEigenvalueField().sage()

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        raise NotImplementedError()

    def trace_of_frobenius(self, prime, power=Integer(1) ):
        """Give the trace of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The trace of $\rho(F_p^n)$, where $\rho$ is the mod $l$ or
        l-adic galois representation associated to this newform, $F_p$
        is the frobenius element at the given prime, and $n$ is the
        given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: K = QuadraticField(2)
            sage: nf = get_newforms(K.ideal(11), base_field=K, algorithm='magma')[0]
            sage: Kf = nf.base_field()
            sage: nf.trace_of_frobenius(Kf.prime_above(7))
            -2
            sage: nf.trace_of_frobenius(Kf.prime_above(17))
            -2
            sage: nf.trace_of_frobenius(Kf.prime_above(23))
            -1
            sage: nf.trace_of_frobenius(Kf.prime_above(31))
            7

        .. SEE_ALSO::

            :meth:`determinant_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        F = self.base_field()
        if not is_Ideal(prime):
            prime = F.ideal(prime)
        if not (prime in F.ideal_monoid()) or not prime.is_prime():
            raise ValueError("%s is not a prime ideal."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == Integer(1) :
            F = self._M.BaseField()
            OF = F.IntegerRing()
            gens = [F(g.list()) for g in prime.gens()]
            prime = None
            for g in gens:
                if prime is None:
                    prime = g*OF
                else:
                    prime += g*OF
            prime = prime.Support().Random()
            return self._f.HeckeEigenvalue(prime).sage()
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=Integer(1) ):
        """Give the determinant of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The determinant of $\rho(F_p^n)$, where $\rho$ is the mod $l$
        or $l$-adic galois representation associated to this newform,
        $F_p$ is the frobenius element at the given prime, and $n$ is
        the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.determinant_of_frobenius(5)
            5
            sage: nf.determinant_of_frobenius(7)
            7
            sage: nf.determinant_of_frobenius(7, power=2)
            49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        raise NotImplementedError()

    def _repr_(self):
        """Give a string representation of this newform"""
        return str(self._f)

    def _latex_(self):
        """Give a latex representation of this newform."""
        return latex(self._f)

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_magma_hilbert(self._M)

class WrappedNewform_magma_bianchi(WrappedNewform):
    r"""A wrapper class around a magma Bianchi modular newform of parallel
    weight 2.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    EXAMPLE::

        sage: from modular_method.modular_forms.newform_wrapper import get_newforms
        sage: K = QuadraticField(-2)
        sage: nf = get_newforms(K.ideal(11), base_field=K, algorithm='magma')[0]; nf
        Element of Cuspidal newform space of Bianchi modular forms
        sage: nf.level()
        Fractional ideal (11)
        sage: nf.base_field()
        Number Field in a with defining polynomial x^2 + 2

    """

    def __init__(self, space):
        r"""Initialize a wrapped newform

        INPUT:

        - ``newform`` -- A Bianchi modular newform as a magma object

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import WrappedNewform_magma_bianchi
            sage: K = QuadraticField(-1)
            sage: cfs = magma.BianchiCuspForms(K, K.ideal(18))
            sage: nfs = cfs.NewSubspace()
            sage: WrappedNewform_magma_bianchi(nfs.NewformDecomposition()[1])
            Element of Cuspidal newform space of Bianchi modular forms

        """
        self._f = space.Eigenform()
        self._M = space
        WrappedNewform.__init__(self)

    def level(self):
        r"""Give the level of this newform.

        OUTPUT:

        An ideal of the totally real number field associated to this
        Hilbert modular form that describes the level of this newform.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.level()
            19

        """
        return self._M.Level().sage()

    @cached_method
    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        return DirichletGroup(Integer(1) )[Integer(0) ]

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        raise NotImplementedError();

    def base_field(self):
        r"""Give the base field for this newform.

        For classical newforms this is always the rationals. For
        Hilbert modular forms this is the totally real field for which
        this is a Hilbert modular form. For Bianchi modular forms it
        is the complex field for which this is a Bianchi modular form.

        OUTPUT:

        The rational field if this newform is a classical modular form.

        The totally real field for which this newform is a Hilbert
        modular form if it is a Hilbert modular form.

        The complex field for which this newform is a Bianchi modular
        form if it is a Bianchi modular form.

        """
        return self._M.BaseField().sage()

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the coefficients of the q-expansion of
        this newform at infinity are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: K = QuadraticField(-3)
            sage: nf = get_newforms(K.ideal(52), base_field=K, algorithm='magma')[1]
            sage: nf.coefficient_field()
            Number Field in K1 with defining polynomial x^2 - 2*x - 4

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._M.HeckeEigenvalueField().sage()

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        raise NotImplementedError()

    def trace_of_frobenius(self, prime, power=Integer(1) ):
        """Give the trace of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The trace of $\rho(F_p^n)$, where $\rho$ is the mod $l$ or
        l-adic galois representation associated to this newform, $F_p$
        is the frobenius element at the given prime, and $n$ is the
        given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: K = QuadraticField(-3)
            sage: nf = get_newforms(K.ideal(26), base_field=K, algorithm='magma')[0]
            sage: Kf = nf.base_field()
            sage: nf.trace_of_frobenius(Kf.prime_above(5))
            -1
            sage: nf.trace_of_frobenius(Kf.prime_above(7))
            -1
            sage: nf.trace_of_frobenius(Kf.prime_above(11))
            14

        .. SEE_ALSO::

            :meth:`determinant_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        F = self.base_field()
        if not is_Ideal(prime):
            prime = F.ideal(prime)
        if not (prime in F.ideal_monoid()) or not prime.is_prime():
            raise ValueError("%s is not a prime ideal."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == Integer(1) :
            F = self._M.BaseField()
            OF = F.IntegerRing()
            gens = [F(g.list()) for g in prime.gens()]
            prime = None
            for g in gens:
                if prime is None:
                    prime = g*OF
                else:
                    prime += g*OF
            prime = prime.Support().Random()
            return self._f.HeckeEigenvalue(prime).sage()
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=Integer(1) ):
        """Give the determinant of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The determinant of $\rho(F_p^n)$, where $\rho$ is the mod $l$
        or $l$-adic galois representation associated to this newform,
        $F_p$ is the frobenius element at the given prime, and $n$ is
        the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.determinant_of_frobenius(5)
            5
            sage: nf.determinant_of_frobenius(7)
            7
            sage: nf.determinant_of_frobenius(7, power=2)
            49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        raise NotImplementedError()

    def _repr_(self):
        """Give a string representation of this newform"""
        return str(self._f)

    def _latex_(self):
        """Give a latex representation of this newform."""
        return latex(self._f)

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_magma_bianchi(self._M)

class WrappedNewform_twisted(WrappedNewform):
    """A wrapper class representing another wrapped newform twisted by a
    character.

    """

    def __init__(self, *, newform, twist):
        r"""Initialize this object"""
        super().__init__()
        self.newform = newform
        self.twist = twist
        L1 = self.twist.base_ring()
        L2 = self.newform.coefficient_field()
        L, phi, psi = common_embedding_field(L1, L2, give_maps=True)
        self._coefficient_field = L
        self._twist_embedding = phi
        self._newform_embedding = psi

    def character(self):
        r"""Give the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform as a
        primitive character.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: eps = DirichletGroup(16).gens()[1]
            sage: nf = get_newforms(16, character=eps)[0]
            sage: nf.character()
            Dirichlet character modulo 16 of conductor 16 mapping 15 |--> 1, 5 |--> zeta4

        """
        return dirichlet_product(self.newform.character(), self.twist)

    def base_field(self):
        r"""Give the base field for this newform.

        For classical newforms this is always the rationals. For
        Hilbert modular forms this is the totally real field for which
        this is a Hilbert modular form. For Bianchi modular forms it
        is the complex field for which this is a Bianchi modular form.

        OUTPUT:

        The rational field if this newform is a classical modular form.

        The totally real field for which this newform is a Hilbert
        modular form if it is a Hilbert modular form.

        The complex field for which this newform is a Bianchi modular
        form if it is a Bianchi modular form.

        """
        return self.newform.base_field()

    def can_compute_frobenius(self, prime):
        r"""Determine whether the trace of Frobenius can be computed for `prime`.

        INPUT:

        - ``prime`` -- A finite prime of the :meth:`base_field`. If
          :meth:`base_field` is the rationals it should be a prime
          number, otherwise a prime ideal.

        OUTPUT:

        - `True`, if this wrapped newform can compute the trace of
          frobenius at the prime `prime`. At least the `prime` should
          then not divide the :meth:`level`.

        - `False`, otherwise.

        """
        return self.newform.can_compute_frobenius(prime) and not prime.divides(self.twist.modulus())

    def coefficient(self, n):
        r"""Give the n-th coefficient of this newform.

        INPUT:

        - ``n`` -- A non-negative integer.

        OUTPUT:

        The n-th coefficient of the q-expansion of this newform at
        infinity.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient(19)
            1
            sage: nf.coefficient(27)
            4

        .. SEE_ALSO::

            :meth:`coefficient_field`,
            :meth:`q_expansion`

        """
        return (self._newform_embedding(self.newform.coefficient(n)) *
                self._twist_embedding(self.twist(n)))

    def coefficient_field(self):
        r"""Give the field over which the coefficients of this newform are
        defined.

        OUTPUT:

        The field over which the eigenvalues of the Hecke operators on
        this newform are defined.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.coefficient_field()
            Rational Field
            sage: nf = get_newforms(31)[0]
            sage: nf.coefficient_field()
            Number Field in a0 with defining polynomial x^2 - x - 1

        .. SEE_ALSO::

            :meth:`coefficient`,
            :meth:`q_expansion`

        """
        return self._coefficient_field

    def trace_of_frobenius(self, prime, power=Integer(1) ):
        """Give the trace of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The trace of $\rho(F_p^n)$, where $\rho$ is the mod $l$ or
        l-adic galois representation associated to this newform, $F_p$
        is the frobenius element at the given prime, and $n$ is the
        given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.trace_of_frobenius(2)
            0
            sage: nf.trace_of_frobenius(7)
            -1
            sage: nf.trace_of_frobenius(7, power=2)
            -13

        .. SEE_ALSO::

            :meth:`determinant_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        if power == Integer(1) :
            p = prime if prime in ZZ else self.base_field.ideal(prime).smallest_integer()
            return (self._newform_embedding(self.newform.trace_of_frobenius(prime)) *
                    self._twist_embedding(self.twist(p)))
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(),
                                          give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=Integer(1) ):
        """Give the determinant of frobenius under the galois representation
        associated to this newform.

        Will give a ValueError if the given prime divides the level of
        this newform, since in that case all mentioned galois
        representations are ramified.

        INPUT:

        - ``prime`` -- A prime of the base field for this newform
          indicating the Frobenius element to be used. This must be a
          prime number if the base field is the rationals and a prime
          ideal otherwise.

        - ``power`` -- A non-negative number (default: 1). If set to
          any value larger than 1, will compute the trace of the
          frobenius element to the given power instead.

        OUTPUT:

        The determinant of $\rho(F_p^n)$, where $\rho$ is the mod $l$
        or $l$-adic galois representation associated to this newform,
        $F_p$ is the frobenius element at the given prime, and $n$ is
        the given argument `power`.

        Since the result does not depend on the choice of $l$, this
        result will be an element of the coefficient field of the
        newform. The only condition is that $l$ and $p$ must be
        coprime.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.determinant_of_frobenius(5)
            5
            sage: nf.determinant_of_frobenius(7)
            7
            sage: nf.determinant_of_frobenius(7, power=2)
            49

        .. SEE_ALSO::

            :meth:`trace_of_frobenius`,
            :meth:`characteristic_polynomial`

        """
        D1 = self.newform.determinant_of_frobenius(prime)
        p = prime if prime in ZZ else self.base_field().ideal(prime).smallest_integer()
        D2 = self.twist(p)
        K, D2_map, D1_map = composite_field(D2.parent(), D1.parent(), give_maps=True)
        D = D1_map(D1) * D2_map(D2)
        return D**power

    def has_cm(self, proof=True):
        """Determine if this newform has complex multiplication.

        INPUT:

        - ``proof`` -- A boolean (default: True). If set to True the
          answer will have been proven correct. If set to False may
          use bounds that have not been proved.

        OUTPUT:

        True if the abelian variety corresponding to this newform has
        complex multiplication. False in any other case.

        EXAMPLE::

            sage: from modular_method.modular_forms.newform_wrapper import get_newforms
            sage: nf = get_newforms(19)[0]
            sage: nf.has_cm()
            False

        """
        return self.newform.has_cm()

    def twist(self, character):
        r"""Twist this newform with a character.

        INPUT:

        - ``character`` -- A Dirichlet character.

        OUTPUT:

        A newform that is this newform twisted by `character`.

        """
        twist = dirichlet_product(self.twist, character)
        if twist.conductor() == Integer(1) :
            return self.newform
        else:
            return WrappedNewform_twisted(newform=self.newform,
                                          twist=twist.primitive_character())

    def copy(self):
        r"""Create a copy of this newform

        The copy saves on memory as much as possible, but is not
        completely shallow as the embeddings list is initialized
        again.

        OUTPUT:

        A copy of this newform.

        """
        return WrappedNewform_twisted(newform=self.newform, twist=self.twist)

    def _repr_(self):
        """Give a string representation of this newform"""
        return (str(self.newform) +
                " twisted by " +
                str(self.twist))

