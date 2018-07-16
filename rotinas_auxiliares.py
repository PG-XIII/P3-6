def prod_interno(vec1, vec2):
    return (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2])

def prod_vetorial(vec1, vec2):
    return [vec1[1]*vec2[2] - vec1[2]*vec2[1], vec1[2]*vec2[0] - vec1[0]*vec2[2], vec1[0]*vec2[1] - vec1[1]*vec2[0]]

def proj_vetores(vec1, vec2):
    prop = prod_interno(vec1,vec2)/prod_interno(vec2,vec2)
    return [prop * vec2[0], prop * vec2[1], prop * vec2[2]]

def normalizar(vec):
    mod = pow(prod_interno(vec,vec),1/2)
    return [vec[0]/mod, vec[1]/mod, vec[2]/mod]

def mul_escalar(vec, k):
    return [vec[0] * k, vec[1] * k, vec[1] * k, vec[2] * k]
