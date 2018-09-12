val_key = {}
for key, val in l.cache.iteritems():
    if val in val_key:
        val_key[val].append(key)
    else:
        val_key[val] = [key]
