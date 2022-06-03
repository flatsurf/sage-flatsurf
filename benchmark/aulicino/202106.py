import pyflatsurf
from flatsurf.geometry.pyflatsurf_conversion import to_pyflatsurf
from pyexactreal import ExactReals
from sage.all import sqrt
from surface_dynamics.all import iet


def build_torus_permutation(n):
    import string

    alphabet = string.ascii_lowercase
    top_list = [alphabet[i] + " " for i in range(n)]
    top_str = ""
    for i in top_list:
        top_str += i
    bottom_list = [alphabet[n - 1] + " "] + [alphabet[i] + " " for i in range(n - 1)]
    bottom_str = ""
    for i in bottom_list:
        bottom_str += i
    return top_str[:-1], bottom_str[:-1]


def target(component, L):
    if component.cylinder():
        # This component is a cylinder. No further decomposition needed.
        return True
    if component.withoutPeriodicTrajectory():
        # This component is minimal. Further decomposition will not produce any cylinders.
        return True

    height = component.height()

    # This height bounds the size of any cylinder. However, it is stretched by the length of the vector
    # defining the vertical direction. (That vector is not normalized because that is hard to do in
    # general rings…)
    bound = (height * height) / pyflatsurf.flatsurf.Bound.upper(
        component.vertical().vertical()
    ).squared()
    return bound > L


def sum_vectors(ray_list):
    sum_vect = []
    for j in enumerate(ray_list[0]):
        sum_1 = sum([i[j[0]] for i in ray_list])
        sum_vect += [sum_1]
    return sum_vect


def construct_cyl_area_list(surface, L):
    def length_period(v):
        return sqrt(float(v.x()) ** 2 + float(v.y()) ** 2)

    directions = S.connections().bound(L).slopes()
    circumferences = []
    areas = []
    # This is twice the surface area, but it washes out in the proportion
    surface_area = float(S.area())
    for direction in directions:
        print(direction)
        from pyflatsurf import flatsurf

        decomposition = flatsurf.makeFlowDecomposition(S, direction.vector())
        decomposition.decompose(lambda component: target(component, L))
        for component in decomposition.components():
            if component.cylinder():
                circumference = component.circumferenceHolonomy()
                if length_period(circumference) > L:
                    continue
                circumferences.append(length_period(circumference))
                # The double area issue disappears here
                areas.append(float(component.area()) / (surface_area))
    return (areas, max(circumferences))


base = iet.Permutation(build_torus_permutation(3))

base_cover = base.cover(["(0, 2, 1)", "()", "(0, 2, 1)"], as_tuple=True)
ray_list_raw = base.suspension_cone().rays_list()
hor_lengths = [ExactReals().random_element() for j in ray_list_raw[0]]
hor_lengths[-1] = ExactReals().one()
S_raw_1 = base_cover.masur_polygon(hor_lengths, sum_vectors(ray_list_raw))
S_raw_2 = to_pyflatsurf(S_raw_1)
S = S_raw_2.eliminateMarkedPoints().surface()

bound = 20
area_list_plus_max = construct_cyl_area_list(S, bound)
