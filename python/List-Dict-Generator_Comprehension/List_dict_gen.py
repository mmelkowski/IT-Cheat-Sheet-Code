# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:17:16 2019

@author: mmelkowski
"""

######################
# List Comprehension #
######################
# ----- # List of pair element
pair = [elt for elt in range(10) if elt % 2 == 0]

# ----- # List of string of pair element
str_pair = [str(elt) for elt in pair]


###########################
# Generator Comprehension #
###########################

# Similaire a la comprehension de liste.
# Permet d'executer la comprehension de list  sans stocker dans la mémoire
# Pratique lorsque le résultat doit être vu une seule fois

# ----- # Generator of pair element
gen_pair = (elt for elt in range(10) if elt % 2 == 0)

# ----- # List of string of pair element
gen_str_pair = (str(elt) for elt in gen_pair)


######################
# Dict Comprehension #
######################
x = {k: v for (k, v) in zip(range(10), range(10))}

d = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
z = {x: d[x] for x in d}
# iteration on dict only give the key

# ----- # Use of condition
y = {x: d[x] for x in d if d[x] > 1}
