from eroML.ensemble import fake_ensemble 

e0 = fake_ensemble(N=11, random_pos=False, center=(10,45), pm=10, seed=1, random_cols=["x","y"])
e1 = fake_ensemble(N=11, random_pos=False, center=(10,45), pm=10, seed=1, random_cols=["a","b"])

#e1 = fake_ensemble(N=11, random_pos=False, center=(1,5), pm=10, seed=1, random_cols=["a","b"])

#print(e0.skyCoords().ra.degree)

e0.merge_add(e1)

e0.append(e1, cols='all', postfix="_NN2")
