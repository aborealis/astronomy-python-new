"""
Astrological calculations for primary directions
"""
from datetime import datetime
from typing import Optional
from components.vector import coords as crd
from sphere import Sphere, true_distance


def possible_aspects(aspect: Optional[int] = None):
    """
    Returns the list of possible aspects (sinister and
    dexter aspects) based on aspect's angle in degrees
    """
    if aspect is None:
        return [-120, -90, - 60, 0, 60, 90, 120, 180]
    if aspect in [60, 90, 120, 180]:
        return [-1 * aspect, aspect]
    return [0]


class Directions:
    """
    Calculates different types of primary directions
    """

    def __init__(self, sphere: Sphere) -> None:
        self.sphere = sphere

    def quadrant(self, right_asc: float, dec: float) -> int:
        """
        Returns the celestial sphere quadrant,
        which a given point belongs to.
        0 - east-top, or east horizon, or MC
        1 - east-bottom
        2 - west-bottom
        3 - west-top, or west horizon
        """
        xyz_hrz = self.sphere.set_equatorial(right_asc, dec).horizontal_xyz()
        if xyz_hrz.x >= 0 and xyz_hrz.z >= 0:
            return 0
        if xyz_hrz.x > 0 and xyz_hrz.z < 0:
            return 1
        if xyz_hrz.x < 0 and xyz_hrz.z < 0:
            return 2
        return 3

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Placidus direction common methods:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def md_portion(self, right_asc: float, dec: float) -> Optional[float]:
        """
        Returns the meridian distance portion,
        MDP of upper/lower meridian distance
        relative to diurnal/nocturnal semiarc
        for a given point
        """
        dsa = self.sphere.dsa(dec)
        umd = self.sphere.umd(right_asc)
        if dsa is None:
            return None

        if umd <= dsa:
            return umd / dsa
        return (180 - umd) / (180 - dsa)

    def __ra_conj_placidus(self,
                           mdp: float,
                           quadrant: int,
                           dec: float) -> float:
        """
        Returns the RA of point A with
        declination dec(A), which `conjuncts`
        point B with the given meridian
        distance portion MDP(B) in a given
        quadrant(B) in placidus system
        """
        ramc = self.sphere.ramc
        dsa = self.sphere.dsa(dec)
        if dsa is None:
            return None
        nsa = 180 - dsa

        if quadrant == 0:
            return ramc + mdp * dsa
        if quadrant == 3:
            return ramc - mdp * dsa
        if quadrant == 1:
            return ramc + 180 - mdp * nsa
        return ramc + 180 + mdp * nsa

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Placidus mundane directions:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def mundane_positions_placidus(self,
                                   acceptor: crd.Vector,
                                   aspect: Optional[int] = None) -> Optional[list[dict]]:
        """
        Returns mundane positions of the acceptor.
        It is the point on equatorial plane with the same
        meridian distance portion in the same quadrant
        as the acceptor, plus all major aspects to
        that point on the equatorial plane.
        """

        acceptor_rasc = acceptor.equatorial().rasc
        acceptor_dec = acceptor.equatorial().dec

        acceptor_quadrant = self.quadrant(acceptor_rasc, acceptor_dec)
        acceptor_mdp = self.md_portion(acceptor_rasc, acceptor_dec)
        if acceptor_mdp is None:
            return None

        # Find RA of the equator plane with the same
        # meridian distance portion as the acceptor
        eqt_rasc = self.__ra_conj_placidus(
            acceptor_mdp, acceptor_quadrant, dec=0
        )

        # Find aspect positions of that point on the equator plane
        return [
            dict(rasc=(eqt_rasc + _a) % 360, aspect=abs(_a))
            for _a in possible_aspects(aspect)
        ]

    def placidus_mundane(self,
                         promissor: crd.Vector,
                         acceptor: crd.Vector,
                         aspect: Optional[int] = None) -> Optional[list[dict]]:
        """
        Returns mundane primary directions of the promissor
        to aspect points of the acceptor
        """
        mundane_positions = self.mundane_positions_placidus(acceptor, aspect)
        if mundane_positions is None:
            return None
        promissor_dec = promissor.equatorial().dec
        promissor_rasc = promissor.equatorial().rasc

        # Find conjunctions of the promissor with these aspect positions
        directions = []
        for item in mundane_positions:
            eqt_quadrant = self.quadrant(item['rasc'], dec=0)
            eqt_mdp = self.md_portion(item['rasc'], dec=0)
            aspect_rasc = self.__ra_conj_placidus(
                eqt_mdp, eqt_quadrant, promissor_dec
            )
            distance = true_distance(aspect_rasc, promissor_rasc)
            check = true_distance(aspect_rasc, promissor_rasc - distance)

            # Negative distances for converse direction
            directions.append(dict(
                dist=distance if abs(check) < 1e-10 else -distance,
                aspect=item['aspect']
            ))

        return directions


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':
    # Initiate a celestial sphere object
    test_sphere = Sphere(
        datetime(2022, 10, 10, 4, 20),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=66 + 43/60,
    )

    test_directions = Directions(test_sphere)

    # Set the 1st point
    proms = test_sphere.set_ecliptical(150, 0)

    # Set the 2nd point
    accpt = test_sphere.set_ecliptical(116.92, 0)

    print(test_directions.placidus_mundane(
        proms, accpt, 0
    ))
