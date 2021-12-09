from typing import List
import hashlib

hash_tree = ['ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb', 
            '3e23e8160039594a33894f6564e1b1348bbd7a0088d42c4acb73eeaed59c009d',
            '2e7d2c03a9507ae265ecf5b5356885a53393a2029d241394997265a1a25aefc6',
            '']

def return_root(hash_tree: List[str])->str:
    hash_tree_2 = []
    i = 0
    while len(hash_tree_2)*2 != len(hash_tree):
        
        st = hash_tree[i] + hash_tree[i + 1]
        i = i + 2
        hash_tree_2.append(hashlib.sha256((st.encode())).hexdigest())
    hash_tree = hash_tree_2.copy()
    if len(hash_tree) == 1:
        return hash_tree[0]
    else:
        return return_root(hash_tree)

