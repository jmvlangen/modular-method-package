def roots_in_field(K):
    r"""
    Gives the elements of $\Z$ that have a root in the given field.

    INPUT:

    - ``K`` -- A number field

    OUTPUT:
    
    A list of all the square free elements of $\Z$ that have
    a root in K.
    """
    return [tmp[0].discriminant().squarefree_part() for tmp in K.subfields(degree=2)]
