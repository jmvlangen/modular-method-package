### Single test
import time

def two_isogeny_test(k):
    start = time.time()
    try:
        result = []
        E = Qcurve_with_2_isogeny(k)
        result.append(E)
        if not E.does_decompose():
            E = E.decomposable_twist()
            result.append(E)
        result.append(E.newform())
    except Exception as e:
        result.append(e)
    end = time.time()
    result.append(end - start)
    return result

### Running tests and saving results:
results = dict()
for i in range(1000):
    results[i] = two_isogeny_test(i)
    print i
    print results[i]
    print ""
    results[-i] = two_isogeny_test(-i)
    print -i
    print results[-i]
    print ""
