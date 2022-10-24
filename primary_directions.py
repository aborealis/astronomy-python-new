"""
Astrological calculations for primary directions
"""
from datetime import datetime
from typing import Optional
from components.vector import coords as crd
from sphere import Sphere, true_distance


class Directions:
    """
    Calculates different types of primary directions
    """

    def __init__(self, sphere: Sphere) -> None:
        self.sphere = sphere
        self.__promissor = None
        self.__acceptor = None

    @property
    def promissor(self):
        """The getter function for promissaor"""
        return self.__promissor

    @property
    def acceptor(self):
        """The getter function for significator"""
        return self.__acceptor

    @promissor.setter
    def promissor(self, coords: crd.Equatorial) -> None:
        if not isinstance(coords, crd.Equatorial):
            raise ValueError(
                "Promissor should be set in equatorial coordinates"
            )
        self.__promissor = coords

    @acceptor.setter
    def acceptor(self, coords: crd.Equatorial) -> None:
        if not isinstance(coords, crd.Equatorial):
            raise ValueError(
                "Acceptor should be set in equatorial coordinates"
            )
        self.__acceptor = coords

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

    def mundane_positions_placidus(self) -> Optional[list[dict]]:
        """
        Returns mundane positions of the acceptor.
        It is the point on equatorial plane with the same
        meridian distance portion in the same quadrant
        as the acceptorm plus all major aspects to
        this point on the equatorial plane.
        """
        if self.__acceptor is None:
            return None

        acceptor_rasc = self.__acceptor.rasc
        acceptor_dec = self.__acceptor.dec

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
            dict(rasc=(eqt_rasc + aspect) % 360, aspect=abs(aspect))
            for aspect in [-120, -90, -60, 0, 60, 90, 120, 180]
        ]

    def direction_placidus(self) -> Optional[list[dict]]:
        """
        Returns primary directions of the promissor
        to aspect points of the acceptor
        """
        mundane_positions = self.mundane_positions_placidus()
        if mundane_positions is None:
            return None
        promissor_dec = self.__promissor.dec
        promissor_rasc = self.__promissor.rasc

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
    point1 = test_sphere.set_ecliptical(150, 0).equatorial()

    # Set the 2nd point
    point2 = test_sphere.set_ecliptical(116.92, 0).equatorial()

    # Init promissor and significator
    test_directions.promissor = point1
    test_directions.acceptor = point2

    print(test_directions.direction_placidus())
