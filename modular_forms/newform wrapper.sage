r"""Tools to work with newforms from different sources.

This file provides a wrapper class that has subclasses wrapping around
different kinds of newforms to provide a uniform way of working with
newforms from different sources.

The base class for wrapped newforms is :class:`Newform_wrapped`. It
has subclasses wrapping around a newform produced by Sage
(:class:`Newform_wrapped_sage`), wrapping around a newform produced by
Magma (:class:`Newform_wrapped_magma`) and wrapping around a newform
produced by reading fourier coefficients from a file
(:class:`Newform_wrapped_stored`).

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

TODO

AUTHORS:

- Joey van Langen (2019-03-01): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Joey van Langen <j.m.van.langen@vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import re

from sage.modular.dirichlet import is_DirichletCharacter

def get_newforms(level, character=None, algorithm='sage', minimal_coeffs=QQ,
                 names='a', path=None):
    r"""Compute the newforms of a given level and character.

    INPUT:

    - ``level`` -- A strictly positive integer which is the level of
      the newforms.

    - ``character`` -- A dirichlet character of which the conductor
      divides the given level.

    - ``algorithm`` -- One of the following possible arguments: 'sage'
      (default) if sage should be used to compute the newforms;
      'magma' if magma should be used to compute the newforms; or
      'file' if the newforms should be loaded from a file.

    - ``minimal_coeffs`` -- A number field or the rationals (default:
      QQ) that should be contained in the coefficient field of each
      newform computed.

    - ``names`` -- An argument required by the sage implementation of
      newforms to be used as the names for the generator of
      coefficient fields of newforms that are not QQ.

    - ``path`` -- A string or None (default: None). Only used in case
      the algorithm is set to file, in which case it should be the
      path to the file from which to load the newforms as a string.
    
    OUTPUT:
    
    A list of instances of Newform_wrapped that contains exactly one
    newform in each galois orbit of newforms in $S(\Gamma_1(N),
    \varepsilon)$, wher $N$ is the given level and $\varepsilon$ is
    the given character. Furthermore the coefficient field of each of
    these newforms extends the given field minimal_coeffs.

    """
    if algorithm == 'sage':
        if character is None:
            nfs = Newforms(level, names=names)
        else:
            eps = character.primitive_character().extend(level)
            nfs = Newforms(eps, names=names)
        result = [Newform_wrapped_sage(f) for f in nfs]
    elif algorithm == 'magma':
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
        result = [Newform_wrapped_magma(orbit[1]) for orbit in nfs]
    elif algorithm == 'file':
        if character is None:
            character = DirichletGroup(1)[0]
        character = character.primitive_character()
        if path is None:
            raise ValueError("Argument path should be set if algorithm file " +
                             "is chosen.")
        to_do = [load_newforms(path)]
        result = []
        while len(to_do) > 0:
            check = to_do
            to_do = []
            for element in check:
                if isinstance(element, list):
                    to_do.extend(element)
                elif (isinstance(element, Newform_wrapped) and
                      element.level() == level and
                      element.character().primitive_character() == character):
                    result.append(element)
                else:
                    pass # Not a right newform, skip!
    else:
        raise ValueError("%s is not a valid algorithm."%(algorithm,))
    if minimal_coeffs == QQ:
        return result
    else:
        poly = minimal_coeffs.gen().minpoly()
        return [f for f in result
                if len(minimal_coeffs.embeddings(f.coefficient_field())) > 0]

class LabeledElement:
    r"""A helper class that stores an element and a label."""
    def __init__(self, label, element):
        self.label = label
        self.element = element
    
def save_newforms(newforms, file_name, coefficients=50, repr_coefficients=True,
                  save_cm=True):
    r"""Save newforms to a file.

    Saves a newform or a list of newforms to a file. This file will
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

    - A newform is represented by a list containing an integer
      preceded by the identifier 'level', a boolean preceded by the
      identifier 'cm', a character, a number field and a list of
      values, all preceded by the identifier 'newform'. The first
      entry will be the level of the newform, the second a boolean
      indicating whether or not this newform has complex multiplication,
      the third the corresponding character, the fourth the coefficient
      field of the newform and the last the coefficients of the newform
      (at some) indices. The entry with label 'cm' may be left out
      or set to -1 to indicate that this information is not known.

    - A list of things will be represented as a list of the
      corresponding representations.

    INPUT:

    - ``newforms`` -- An instance of Newform_wrapped. This may also be
      a list or other iterable containing as elements instances of
      Newform_wrapped or lists that satisfy the same property. These
      are the newforms that will be saved to the file.

    - ``file_name`` -- A string containing the file name to which the
      given newforms should be saved.

    - ``coefficients`` -- An iterable object of non-negative integers
      or a non-negative integer (default: 50). This determines the
      coefficients of the newform that will be saved. If it is an
      iterable object, will save all the coefficients with the indices
      given in that object. If it is a single non-negative integer
      will make this to be the list of all non-negative integers up to
      the given integer (excluding the integer itself).

    - ``repr_coefficients`` -- A boolean value (default: True).  If
      set to true, will always save the first 20 coeficients of this
      newform. This is recommmended as otherwise the string
      representation of the newforms will break if loaded back in.

    - ``save_cm`` -- A boolean value (default: True) indicating
      whether for each newform saved the information whether or not it
      has complex multiplication should also be computed and saved. If
      set to False, this will not be done and the field 'cm' of a
      newform will be set to -1.

    EXAMPLE::

    TODO

    """
    if coefficients in ZZ and coefficients > 0:
        coefficients = range(coefficients)
    if repr_coefficients:
        coefficients = range(20) + [c for c in coefficients if c > 20]
    else:
        coefficients = [c for c in coefficients]
    with open(file_name, "w+") as f:
        _write_element(newforms, f, coefficients, save_cm)

def _write_list(ls, f, coefficients, save_cm, indent=0, indent_start=True):
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
        f.write(" "*4*indent)
    f.write('[')
    write_comma = False
    for element in ls:
        if write_comma:
            f.write(',')
        f.write('\n')
        _write_element(element, f, coefficients, save_cm, indent=indent+1)
        write_comma=True
    f.write(']')

def _write_element(element, f, coefficients, save_cm, indent=0):
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
    if isinstance(element, Newform_wrapped):
        _write_newform(element, f, coefficients, save_cm, indent=indent)
    elif is_DirichletCharacter(element):
        _write_character(element, f, coefficients, save_cm, indent=indent)
    elif is_NumberField(element):
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

def _write_newform(nf, f, coefficients, save_cm, indent=0):
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
    level = LabeledElement('level', nf.level())
    if save_cm:
        cm = LabeledElement('cm', ZZ(nf.has_cm()))
    else:
        cm = LabeledElement('cm', ZZ(-1))
    character = nf.character()
    field = nf.coefficient_field()
    if not field.is_absolute():
        field = field.absolute_field(names='a')
    values = [(n, (QQ(nf.coefficient(n)) if nf.coefficient(n) in QQ
                   else LabeledElement('element',
                                       field(nf.coefficient(n)).list())))
              for n in coefficients]
    values = LabeledElement('values', values)
    element = LabeledElement('newform', [level, cm, character, field, values])
    _write_labeled_element(element, f, coefficients, save_cm, indent=indent)

def _write_character(eps, f, coefficients, save_cm, indent=0):
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

def _write_field(field, f, coefficients, save_cm, indent=0):
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

def _write_polynomial(poly, f, coefficients, save_cm, indent=0):
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

def _write_rational(q, f, indent=0, indent_start=True):
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
        f.write(" "*4*indent)
    f.write(str(QQ(q)))

def _write_labeled_element(element, f, coefficients, save_cm, indent=0):
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
    f.write(" "*4*indent)
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
    r"""
    Loads newforms from a file.

    A file that stores newforms will adhere to the following
    formatting rules in which we shall use:
 
    <list> := '[' ( <element> ( ',' <element> )* )? ']'
    <element> := ( '<' <identifier> '>' ':=' )? ( <list> | <rational> )
    <identifier> := <letter>+
    <rational> := <integer> ( '/' <positive_integer> )?
    <integer> := <zero> | ( ( '-' )?<positive_integer> )
    <positive_integer> := <non_zero_digit> ( <digit> )*
    <digit> := <zero> | <non_zero_digit>
    <non_zero_digit> := [1-9]
    <zero> := '0'
    <letter> := [a-zA-Z]

    - A boolean value is represented by the integer 1 if it is
      True and the integer 0 if it is False.
    - An element of a number field is represented as the list of
      rational coefficients with respect to the power basis in the
      generator, preceded by the identifier 'element'
    - A polynomial with rational coefficients is represented by
      a list of its coefficients (starting at the constant term)
      preceded by the identifier 'polynomial'
    - A number field is represented by a list containing a polynomial,
      preceded by the identifier 'field'. The polynomial is the
      defining polynomial of the number field.
    - A list of values of a function is represented by a list containing
      lists of exactly two elements, all preceded by the identifier
      'values'. The function maps the first element of a list in the
      corresponding list to the second element thereof.
    - A character is represented by a list containing an integer
      preceded by the identifier 'conductor' and a list of values of
      the character, all preceded by the identifier 'character'.
      The integer with identifier 'conductor' will be the conductor
      of the character. The entry labeled values will be pairs of
      integers, such that if zeta is the relevant n-th root of unity
      a pair (k, e) appears in this list if the character takes the
      value zeta^e at k.
    - A newform is represented by a list containing an integer
      preceded by the identifier 'level', a boolean preceded by the
      identifier 'cm', a character, a number field and a list of
      values, all preceded by the identifier 'newform'. The first
      entry will be the level of the newform, the second a boolean
      indicating whether or not this newform has complex multiplication,
      the third the corresponding character, the fourth the coefficient
      field of the newform and the last the coefficients of the newform
      (at some) indices. The entry with label 'cm' may be left out
      or set to -1 to indicate that this information is not known.
    - A list of things will be represented as a list of the
      corresponding representations.

    INPUT:

    - ``file_name`` -- A string containing the file name from which
      the given newforms should be loaded.

    OUTPUT:

    An instance of Newform_wrapped_stored or a list thereof
    representing the newforms found in the given file.
    """
    with open(file_name, 'r') as f:
        return _read_element(f)

def _interpret_element(element):
    if element.label == 'newform':
        return _interpret_newform(element.element)
    elif element.label == 'character':
        return _interpret_character(element.element)
    elif element.label == 'field':
        return _interpret_field(element.element)
    elif element.label == 'polynomial':
        return _interpret_polynomial(element.element)
    else:
        return element

def _interpret_polynomial(element):
    R.<x> = QQ[]
    return R(element)

def _interpret_field(element):
    if not isinstance(element, list):
        raise ValueError("%s is not a list."%(element,))
    if len(element) != 1:
        raise ValueError("%s does not have length 1."%(element,))
    polynomial = element[0]
    if not is_Polynomial(polynomial):
        raise ValueError("%s is not a polynomial."%(polynomial,))
    if polynomial.degree() == 1:
        return QQ
    return NumberField(polynomial, names='a')

def _interpret_character(element):
    if not isinstance(element, list):
        raise valueError("%s is not a list."%(element,))
    conductor=None
    values=None
    for part in element:
        if isinstance(part, LabeledElement):
            if part.label.lower() == 'conductor' and conductor is None and part.element in ZZ:
                conductor = ZZ(part.element)
            elif part.label.lower() == 'values' and values is None and isinstance(part.element, list):
                values=dict()
                for pair in part.element:
                    if not isinstance(pair, list) or len(pair) != 2 or pair[0] not in ZZ:
                        raise ValueError("Expected a pair for character values, but got %s"%(pair,))
                    if pair[1] in ZZ:
                        values[ZZ(pair[0])] = ZZ(pair[1])
                    else:
                        raise ValueError("Expected a pair for character values, but got %s"%(pair,))
            else:
                raise ValueError("Unexpected element %s with label %s for character."%(part.element, part.label))
        else:
            raise ValueError("Unexpected element %s for character."%(part,))
    if conductor is None or values is None:
        raise ValueError("Not enough arguments to make a character.")
    D = DirichletGroup(conductor)
    try:
        pows = D._zeta_powers
        return D([pows[values[n]] for n in D.unit_gens()])
    except KeyError as e:
        raise ValueError("Requires value at %s to construct Dirichlet character."%(str(e),))

def _interpret_newform(element):
    if not isinstance(element, list):
        raise valueError("%s is not a list."%(element,))
    level=None
    cm=None
    character=None
    field=None
    coefficients=None
    for part in element:
        if is_DirichletCharacter(part) and character is None:
            character=part
        elif is_NumberField(part) and field is None:
            field=part
        elif isinstance(part, LabeledElement):
            if part.label.lower() == 'level' and level is None and part.element in ZZ:
                level = ZZ(part.element)
            elif part.label.lower() == 'values' and coefficients is None and isinstance(part.element, list):
                coefficients=dict()
                for pair in part.element:
                    if not isinstance(pair, list) or len(pair) != 2 or pair[0] not in ZZ:
                        raise ValueError("Expected a pair for newform coefficients, but got %s"%(pair,))
                    if pair[1] in QQ:
                        coefficients[ZZ(pair[0])] = QQ(pair[1])
                    elif isinstance(pair[1], LabeledElement) and pair[1].label.lower() == 'element':
                        coefficients[ZZ(pair[0])] = pair[1].element
                    else:
                        raise ValueError("Expected a pair for newform coefficients, but got %s"%(pair,))
            elif part.label.lower() == 'cm' and cm is None and part.element in ZZ:
                if cm != -1:
                    cm = (part.element != 0)
            else:
                raise ValueError("Unexpected element %s with label %s for newform."%(part.element, part.label))
        else:
            raise ValueError("Unexpected element %s for newform."%(part,))
    if level is None or character is None or field is None or coefficients is None:
        raise ValueError("Not enough arguments to make a newform.")
    for key in coefficients:
        coefficients[key] = field(coefficients[key])
    return Newform_wrapped_stored(level, character, field, coefficients, cm)

def _read_element(f):
    whitespace = re.compile('\s')
    label = None
    element = None
    s = f.read(1)
    while len(s) > 0:
        if whitespace.match(s):
            pass
        elif s == '<':
            f.seek(-1,1)
            label = _read_identifier(f)
            _read_colon_is(f)
        elif s == '[':
            f.seek(-1,1)
            element = _read_list(f)
            break
        elif s.isdigit() or s == '-':
            f.seek(-1,1)
            element = _read_rational(f)
            break
        else:
            raise ValueError("Read unexpected character %s while reading an element."%(s,))
        s = f.read(1)
    else:
        raise ValueError("File ended before reading an element.")
    if label is None:
        return element
    else:
        element = LabeledElement(label, element)
        return _interpret_element(element)

def _read_colon_is(f):
    whitespace = re.compile('\s')
    s = f.read(1)
    while len(s) > 0:
        if whitespace.match(s):
            pass
        elif s == ':':
           s = f.read(1)
           if s == '=':
               return
           else:
               raise ValueError("Attempting to read ':=' at :%s, but failed."%(s,))
        else:
           raise ValueError("Attempting to read ':=' at %s, but failed."%(s,))
        s = f.read(1)
    raise ValueError("Attempting to read ':=', but encountered end of file.")

def _read_identifier(f):
    s = f.read(1)
    if s != '<':
        raise ValueError("Attempting to read '<', but read %s"%(s,))
    result = ""
    s = f.read(1)
    while len(s) > 0:
        if s.isalpha():
            result = result + s
        elif s == '>':
            return result
        else:
            raise ValueError("Attempting to read a letter, but read %s"%(s,))
        s = f.read(1)
    raise ValueError("Reached end of file whilst reading an identifier.")

def _read_list(f):
    whitespace = re.compile('\s')
    s = f.read(1)
    if s != '[':
        raise ValueError("Attempting to read the start of a list at %s, but failed"%(s,))
    result = []
    s = f.read(1)
    can_close = True
    need_comma = False
    while len(s) > 0:
        if whitespace.match(s):
            pass
        elif s == ']' and can_close:
            return result
        elif s == ',' and need_comma:
            need_comma = False
            can_close = False
        elif not need_comma:
            f.seek(-1,1)
            result.append(_read_element(f))
            need_comma = True
            can_close = True
        else:
            raise ValueError("Attempting to read more of a list, but encountered %s"%(s,))
        s = f.read(1)
    raise ValueError("Encountered end of file whilst reading a list.")

def _read_rational(f):
    result = _read_integer(f)
    s = f.read(1)
    if s == '/':
        result = result + '/' + _read_positive_integer(f)
    else:
        f.seek(-1,1)
    return QQ(result)

def _read_integer(f):
    result = ''
    s = f.read(1)
    if s == '-':
        result = result + '-'
    else:
        f.seek(-1,1)
    try:
        result = result + _read_zero(f)
    except ValueError:
        f.seek(-1,1)
        result = result + _read_positive_integer(f)
    return result

def _read_positive_integer(f):
    result = _read_non_zero_digit(f)
    try:
        while True:
            result = result + _read_digit(f)
    except ValueError:
        f = f.seek(-1,1)
        return result
        
def _read_digit(f):
    s = f.read(1)
    if not s.isdigit():
        raise ValueError("Attempting to read a digit, but read %s"%(s,))
    return s
    
def _read_non_zero_digit(f):
    s = f.read(1)
    if s.isdigit() and s != 0:
        return s
    else:
        raise ValueError("Attempting to read a non-zero digit, but read %s"%(s,))
    
def _read_zero(f):
    s = f.read(1)
    if s != '0':
        raise ValueError("Attempting to read 0, but read %s"%(s,))
    return s
    
class Newform_wrapped(SageObject):
    r"""
    The wrapped version of a newform.

    This acts as a common interface to work with a newform,
    independent of its internal representation.

    This class is a base class, but should not be used
    as an instance. It rather provides a template for all
    classes that inherit it.
    """

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        raise NotImplementedError()

    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        raise NotImplementedError()
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        raise NotImplementedError()

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        raise NotImplementedError()

    def q_expansion(self, prec=20):
        """
        Gives the q-expansion of this newform.

        INPUT:
        
        - ``prec`` -- A non-negative integer (default: 20)
          giving a bound on the powers that should be present
          in this q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given
        as a power series in q with coefficients in the
        coefficient field of this newform and capped at
        the prec.
        """
        R.<q> = self.coefficient_field()[[]]
        result = sum(self.coefficient(n) * q^n for n in range(prec))
        return result.add_bigoh(prec)

    @cached_method
    def _trace_power_formula(self, power):
        R.<x,y> = QQ[]
        return polynomial_to_symmetric(x^power + y^power)

    def trace_of_frobenius(self, prime, power=1):
        """
        Gives the trace of frobenius under the galois
        representation associated to this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The trace of rho_f(frob_p^n), where rho_f is the
        mod-l or l-adic galois representation associated
        to this newform, frob_p is the frobenius element
        at prime, and n is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        if power == 1:
            return self.coefficient(prime)
        T = self.trace_of_frobenius(prime)
        D = self.determinant_of_frobenius(prime)
        K, T_map, D_map = composite_field(T.parent(), D.parent(), give_maps=True)
        return self._trace_power_formula(power)(T_map(T), D_map(D))

    def determinant_of_frobenius(self, prime, power=1):
        """
        Gives the determinant of frobenius under the
        galois representation associated to this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The determinant of rho_f(frob_p^n), where rho_f is the
        mod-l or l-adic galois representation associated
        to this newform, frob_p is the frobenius element
        at prime, and n is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        if prime not in ZZ or not prime.is_prime():
            raise ValueError("%s is not a prime number."%(prime,))
        if prime.divides(self.level()):
            raise ValueError("%s divides the level: %s."%(prime, self.level()))
        D = self.character()(prime) * prime
        return D^power

    def characteristic_polynomial(self, prime, power=1):
        """
        Gives the characteristic polynomial of the frobenius
        element acting on this newform.
        
        INPUT:

        - ``prime`` -- A prime number indicating the
          frobenius element to be used.
        - ``power`` -- A non-negative number (default: 1)
          If set to any value larger than 1, will compute
          the trace of the frobenius element to the given
          power instead.

        OUTPUT:

        The characteristic polynomial of rho_f(frob_p^n),
        where rho_f is the mod-l or l-adic galois
        representation associated to this newform,
        frob_p is the frobenius element at prime, and n
        is the given argument power.

        Since the result does not depend on the choice
        of l, this result will be an element of the
        coefficient field of the newform. The only
        condition is that l and p must be coprime.

        Will give a ValueError if the prime divides the
        level of this newform, since in that case all
        mentioned galois representations are ramified.
        """
        T = self.trace_of_frobenius(prime, power=power)
        D = self.determinant_of_frobenius(prime, power=power)
        K, T_map, D_map = composite_field(T.parent(), D.parent(), give_maps=True)
        R.<x> = K[]
        return x^2 - T_map(T)*x + D_map(D)

    def has_cm(self):
        """
        Determines if this newform has complex
        multiplication.

        OUTPUT:

        True if the abelian variety corresponding
        to this newform has complex multiplication.
        False in any other case.
        """
        raise NotImplementedError()
    
    qexp = q_expansion

    def _repr_(self):
        """
        Gives a string representation of self
        """
        return str(self.q_expansion())

    def _latex_(self):
        """
        Gives a latex representation of self.
        """
        return latex(self.q_expansion())

class Newform_wrapped_sage(Newform_wrapped):
    r"""
    The wrapped version of a newform in sage.
    """

    def __init__(self, newform):
        self._f = newform

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        return self._f.level()

    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        return self._f.character()
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        if n == 0:
            return self.coefficient_field()(0)
        else:
            return self._f.coefficient(n)

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        return self._f.base_ring()

    def q_expansion(self, prec=20):
        """
        Gives the q-expansion of this newform.

        INPUT:
        
        - ``prec`` -- A non-negative integer (default: 20)
          giving a bound on the powers that should be present
          in this q-expansion.

        OUTPUT:

        The q-expansion of this newform at infinity given
        as a power series in q with coefficients in the
        coefficient field of this newform and capped at
        the prec.
        """
        return self._f.q_expansion(prec=20)

    def has_cm(self):
        """
        Determines if this newform has complex
        multiplication.

        OUTPUT:

        True if the abelian variety corresponding
        to this newform has complex multiplication.
        False in any other case.
        """
        return self._f.has_cm()

    def _repr_(self):
        """
        Gives a string representation of self.
        """
        return str(self._f)

    def _latex_(self):
        """Gives a latex representation of self."""
        return latex(self._f)

    qexp = q_expansion

class Newform_wrapped_magma(Newform_wrapped):
    r"""
    The wrapped version of a newform in magma.
    """

    def __init__(self, newform):
        self._f = newform

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        return self._f.Level().sage()

    @cached_method
    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
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
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        return self._f.Coefficient(n).sage()

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        return self._f.BaseField().sage()

    def has_cm(self):
        """
        Determines if this newform has complex
        multiplication.

        OUTPUT:

        True if the abelian variety corresponding
        to this newform has complex multiplication.
        False in any other case.
        """
        return self._f.ModularSymbols().HasCM()

    def _repr_(self):
        """Gives a string representation of self."""
        return str(self._f)

    def _latex_(self):
        """Gives a latex representation of self."""
        return latex(self._f)

class Newform_wrapped_stored(Newform_wrapped):
    r"""
    The wrapped version of a newform created from some stored data.
    """

    def __init__(self, level, character, coefficient_field, coefficients, cm):
        self._level = level
        self._eps = character
        self._K = coefficient_field
        self._coeffs = coefficients
        self._cm = cm

    def level(self):
        r"""
        Gives the level of the newform.
        
        OUTPUT:

        A non-negative integer describing the level of this
        newform.
        """
        return self._level

    @cached_method
    def character(self):
        r"""
        Gives the character associated to this newform.

        OUTPUT:

        The dirichlet character associated to this newform
        as a primitive character.
        """
        return self._eps
        
    def coefficient(self, n):
        r"""
        Give the n-th coefficient of this newform.

        INPUT:
        
        - ``n`` -- A non-negative integer.

        OUTPUT:
        
        The n-th coefficient of the q-expansion of
        this newform at infinity.
        """
        try:
            return self._coeffs[n]
        except KeyError:
            raise ValueError("The %s-th coefficient is not stored."%(n,))

    def coefficient_field(self):
        r"""
        Gives the field over which the coefficients of
        this newform are defined.

        OUTPUT:

        The field over which the coefficients of the
        q-expansion of this newform at infinity are
        defined.
        """
        return self._K

    def has_cm(self):
        """
        Determines if this newform has complex
        multiplication.

        OUTPUT:

        True if the abelian variety corresponding
        to this newform has complex multiplication.
        False in any other case.
        """
        if self._cm is None:
            raise ValueError("Undetermined whether this newform has CM.")
        return self._cm

    def _repr_(self):
        """
        Gives a string representation of self
        """
        try:
            return str(self.q_expansion())
        except ValueError:
            return "Loaded newform with limited coefficients."

    def _latex_(self):
        """
        Gives a latex representation of self.
        """
        try:
            return latex(self.q_expansion())
        except ValueError:
            return "\\text{Loaded newform with limited coefficients}"
