# A list of things that would be useful to change/add sorted by module

## modular_method.modular_forms.elimination

* Put apply_to_conditional_value in outermost function everywhere

* Support conditional values in all inputs

  - This would include inside the tuples that can be passed as first
    or second argument

  - Probably requires additional conditional value helpers first, at
    least a simplify on a conditional_product

* Use positional and keyword only arguments

* Replace integer at the end with an object that can register more information

  - In particular it should register primes being eliminated when also
    an infinite number of primes remain, i.e. eliminating l even
    though B = 0

* Significantly speed up the conditions for Kraus method

## modular_method.diophantine_equations.conditions

* Make apply_to_conditional_function a method of ConditionalValue

* Add a modulus to every condition

  - modulus N indicates this condition only depends on the value of
    the parameters modulo Z/N (Z may be ring of integers here)

  - can be used to see for which primes taking the pAdic_tree actually
    is useful

* Add a simplify functions to conditional_values

  - Way to group conditions in long and/or statements based on modulus

  - Way to convert into TreeConditions for non-zero modulus

## modular_method.modular_forms.newform_wrapper

* Method to compute the elliptic curve corresponding to a rational newform

  - supported by Sage/Magma implementation of classical modular forms

## modular_method.elliptic_curves.frey_curves

* Method specialize should not have a single argument that is a tuple

  - use *args instead