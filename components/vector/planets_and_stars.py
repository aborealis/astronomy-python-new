"""
This module calculates the positions of the stars and planets
for a given moment of time
"""
import os
from typing import NamedTuple
import swisseph as swe
from sphere import Sphere


curr_path = os.path.abspath(os.path.dirname(__file__))
swe.set_ephe_path(curr_path)

ALL_STARS = [
    "Abhijit", "Aboras", "Abyssus Aqueus", "Acamar", "Achernar", "Achird", "Acrab", "Acrux",
    "Acubens", "Aculeus", "Acumen", "Adara", "Adhab", "Adhafera", "Adhara", "Adhil", "Agena",
    "Ahadi", "Ain", "Ain al Rami", "Ainalhai", "Akrab", "Al Anz", "Al Athfar", "Al Dhanab",
    "Al Durajah", "Al Fahd", "Al Fawaris", "Al Hamalain", "Al Haud", "Al Hecka", "Al Jabhah",
    "Al Jathiyah", "Al Kalb al Rai", "Al Khabdhilinan", "Al Krikab", "Al Kurud",
    "Al Minliar al Asad", "Al Minliar al Shuja", "Al Mizan", "Al Nitham", "Al Pherg",
    "Al Rescha", "Al Rihla", "Al Sadr al Ketus", "Al Tarf", "Al Thalimaim Anterior",
    "Al Thalimaim Posterior", "Aladfar", "Alagemin", "Alaraph", "Alathfar", "Albaldah",
    "Albali", "Albireo", "Albulaan", "Alchiba", "Alchita", "Alcor", "Alcyone", "Aldafirah",
    "Aldebaran", "Alderamin", "Aldhibah", "Aldhibain", "Alfecca Meridiana", "Alfirk", "Algedi",
    "Algenib", "Algieba", "Algol", "Algorab", "Alhakim", "Alhena", "Alherem",
    "Alifa Al Farkadain", "Alioth", "Alkaid", "Alkalurops", "Alkes", "Alkibash", "Alkidr",
    "Alkurhah", "Almaak", "Almac", "Almach", "Almak", "Almeisan", "Alminhar", "Alnair",
    "Alnasl", "Alnath", "Alnilam", "Alnitak", "Alniyat", "Alphard", "Alphecca", "Alphekka",
    "Alpheratz", "Alphirk", "Alradif", "Alrai", "Alredif", "Alrischa", "Alsafi", "Alsciaukat",
    "Alshain", "Alsharasif", "Alshat", "Alsuhail", "Altager", "Altair", "Altais", "Altawk",
    "Alterf", "Althaur", "Aludra", "Alula Australis", "Alula Borealis", "Alvahet", "Alvashak",
    "Alwaid", "Alya", "Alzirr", "Anazitisi", "Ancha", "Andromeda Galaxy", "Angetenar", "Ankaa",
    "Anser", "Antares", "Anunitum", "Anuradha", "Anwar al Farkadain", "Apex", "Ara",
    "Arcturus", "Arkab Posterior", "Arkab Prior", "Armus", "Arneb", "Arrakis", "Ascella",
    "Asellus Australis", "Asellus Borealis", "Asellus Primus", "Asellus Secundus",
    "Asellus Tertius", "Ashlesha", "Ashvini", "Aspidiske", "Asterion", "Asterope", "Athafi",
    "Atik", "Atiks", "Atirsagne", "Atlas", "Atria", "Auva", "Avior", "Avis Satyra",
    "Azelfafage", "Azha", "Azmidiske", "Baham", "Barnard's star", "Baten Algiedi",
    "Baten Kaitos", "Batentaban Australis", "Batentaban Borealis", "Bazak", "Beid",
    "Bellatrix", "Benetnash", "Bered", "Beteigeuse", "Betelgeuse", "Bharani", "Biham",
    "Birdun", "Bogardus", "Bos", "Botein", "Brachium", "Bunda", "Bungula", "Canopus",
    "Capella", "Caph", "Capulus", "Castor", "Castra", "Cauda Hydrae", "Ceginus", "Celbalrai",
    "Celeano", "Cestan", "Chara", "Chertan", "Chin", "Chort", "Chow", "Cih", "Citra",
    "Cor Caroli", "Cor Hydrae", "Cor Serpentis", "Coxa", "Cujam", "Cursa", "Dabih", "Decapoda",
    "Decrux", "Deli", "Deneb", "Deneb Adige", "Deneb Algedi", "Deneb Algenubi", "Deneb Dulphim",
    "Deneb Kaitos", "Deneb el Okab Australis", "Deneb el Okab Borealis", "Denebola",
    "Dhanab al Shuja", "Dhanishtha", "Dheneb", "Dhur", "Diadem", "Difda", "Diphda", "Dorsum",
    "Drus", "Drys", "Dschubba", "Dubhe", "Dziban", "Ed Asich", "Edasich", "Ekkhysis",
    "El Kophrah", "El Nath", "Electra", "Elnath", "Eltanin", "Enif", "Erakis", "Errai",
    "Etamin", "Facies", "Farkadain", "Fomalhaut", "Foramen", "Fornacis", "Fudail",
    "Fum Alsamakah", "Furibundus", "Furud", "GCRS00", "Gacrux", "Gal. Center", "Gal. Center",
    "Gemma", "Ghusn al Zaitun", "Gianfar", "Giansar", "Giedi Prima", "Giedi Secunda",
    "Gienah Corvi", "Gienah Cygni", "Gienah Ghurab", "Girtab", "Gomeisa", "Gorgona Quatra",
    "Gorgona Secunda", "Gorgona Tertia", "Graffias", "Grafias", "Great Attractor", "Gruid",
    "Grumium", "Hadar", "Haedi", "Haedus", "Hamal", "Han", "Haris", "Hasseleh", "Hasta",
    "Hatsya", "Hecatebolus", "Heka", "Helkath", "Hemelein Prima", "Hemelein Secunda", "Heze",
    "Hilasmus", "Hoedus I", "Hoedus II", "Homam", "Hyadum I", "Hyadum II", "Hydor", "Hydria",
    "Hydrobius", "Imad", "Ira Furoris", "Isis", "Izar", "Jabbah", "Jabhat al Akrab",
    "Jabhat al Akrab", "Jih", "Juxta Crucem", "Jyeshtha", "Kabkent Secunda", "Kabkent Tertia",
    "Kaffaljidhma", "Kaht", "Kaimana", "Kajam epsHer", "Kajam omeHer", "Kakkab", "Kattupothu",
    "Kaus Australis", "Kaus Borealis", "Kaus Medis", "Kaus Meridionalis", "Ke Kwan", "Kebash",
    "Keid", "Kekouan", "Kelb Alrai", "Kerb", "Ketu", "Khambalia", "Kissin", "Kitalpha",
    "Kochab", "Koo She", "Kornephoros", "Kraz", "Krttika", "Kuma", "Kuma", "Kurdah", "Kurhah",
    "Labrum", "Leiolepidotus", "Leiolepis", "Lesath", "Linteum", "Ma Ti", "Maasym", "Maaz",
    "Mabsuthat", "Mabsuthat", "Maculata", "Maculosa", "Magha", "Maia", "Manubrium", "Manus",
    "Marakk", "Marfak", "Marfik", "Markab", "Markeb", "Marsik", "Masym", "Matar", "Mati",
    "Mautinah", "Mebsuta", "Megrez", "Meissa", "Mekbuda", "Melkarth", "Men", "Menkalinan",
    "Menkar", "Menkent", "Menkib", "Merak", "Merga", "Merope", "Mesarthim", "Metallah",
    "Miaplacidus", "Mimosa", "Minchir", "Minkar", "Mintaka", "Mira", "Mirach", "Mirak",
    "Miram", "Mirfak", "Mirphak", "Mirzam", "Misam", "Mizar", "Mrgashirsha", "Mufrid",
    "Muhlifain", "Mula", "Muliphein", "Muphrid", "Murzim", "Murzims", "Muscida", "Nageba",
    "Nair al Saif", "Nanto", "Naos", "Nash", "Nashira", "Nehushtan", "Nekkar", "Nihal",
    "Nodus I", "Nodus II", "Nulla Pambu", "Nunki", "Nusakan", "Oculus", "Pazan", "Pazhan",
    "Peacock", "Peregrini", "Phact", "Phaeo", "Phaesula", "Pharkadain", "Phecda", "Pherkad",
    "Pherkad Minor", "Phicares", "Phicareus", "Phycochroma", "Pleione", "Pleura", "Polaris",
    "Polaris Australis", "Polis", "Pollux", "Porrima", "Praecipua", "Praesepe Cluster",
    "Prijipati", "Prima Hyadum", "Princeps", "Procyon", "Propus etaGem", "Propus iotGem",
    "Proxima Centauri", "Pulcherrima", "Punarvasu", "Purvabhadra", "Purvaphalguni",
    "Purvashadha", "Pushya", "Pushya", "Qin", "Rana", "Ras Algethi", "Ras Alkurki",
    "Ras Elased Australis", "Ras Elased Borealis", "Ras Mutallah", "Rasalas", "Rasalgethi",
    "Rasalhague", "Rastaban", "Regor", "Regulus", "Rehla", "Revati", "Rigel", "Rigel Kentaurus",
    "Rigil Kent", "Rijl al Awwa", "Ril Alauva", "Rohini", "Rotanev", "Ruc", "Rucha", "Ruchbah",
    "Ruchbah I", "Ruchbah II", "Rukbalgethi Genubi", "Rukbalgethi Shemali", "Rukbat", "Rukh",
    "Rutilicus", "Sabik", "Sadalachbia", "Sadalbari", "Sadalmelek", "Sadalmelik", "Sadalsuud",
    "Sadaltager", "Sadatoni", "Sadir", "Sador", "Sadr", "Saidak", "Saiph", "Salm", "Samakah",
    "Sanduleak", "Sargas", "Sarin", "Sasin", "Sataghni", "Sceptrum", "Scheat", "Schedar",
    "Schedir", "Scutulum", "Seat", "Secunda Hyadum", "Segin", "Seginus", "Sephdar", "Sham",
    "Shatabhishaj", "Shaula", "Shedir", "Sheliak", "Shemali", "Sheratan", "Shir", "Shishimai",
    "Shravana", "Shravishtha", "Simiram", "Simmah", "Sinistra", "Sirius", "Sirrah", "Situla",
    "Skat", "Sofian", "Spica", "Spica", "Spiculum", "Sterope I", "Sterope II", "Sualocin",
    "Subra", "Suhail", "Suhail Hadar", "Suhail al Muhlif", "Sulafat", "Sulaphat", "Svati",
    "Syrma", "Tabit", "Talitha Australis", "Talitha Borealis", "Tang", "Tania Australis",
    "Tania Borealis", "Tarazed", "Taygeta", "Tegmen", "Tegmine", "Tejat", "Terebellium",
    "Thabit", "The Blaze Star", "The Garnet Star", "Theemin", "Thuban", "Thusia", "Tien Kang",
    "Toliman", "Torcularis Septentrionalis", "Trapezium", "Tse Tseng", "Tseen Foo", "Tseen Ke",
    "Tsih", "Tsze", "Tsze Tseang", "Turais", "Tureis", "Tyl", "Ukdah", "Unukalhai", "Urakhga",
    "Urodelus", "Ushakaron", "Uttarabhadra", "Uttaraphalguni", "Uttarashadha",
    "Vathorz Posterior", "Vathorz Prior", "Vega", "Vernalis", "Vindemiatrix", "Virgo Cluster",
    "Vishakha", "Wasat", "Wazn", "Wei", "Wezen", "Xestus", "Yed Posterior", "Yed Prior",
    "Yildun", "Zaniah", "Zaurak", "Zavijava", "Zhou", "Zibal", "Zosma", "Zuben Elakrab",
    "Zuben Elakribi", "Zuben Elgenubi", "Zuben Eschamali", "Zuben Hakrabi", "Zubenelakrab",
    "Zubenelakribi",
]


class Planet(NamedTuple):
    """
    Planetary zodiacal coordinates
    """
    lon: float
    lat: float


class Star(NamedTuple):
    """
    A fixes star zodiacal coordinates
    """
    name: str
    lon: float
    lat: float


class StellarObject:
    """
    Contains methods to calculate positions
    of the stars and planets on a given
    celestial sphere.
    """

    def __init__(self, sphere: Sphere) -> None:
        self.jday = sphere.jday
        self.planets = self.__scan_planets__()
        self.stars = self.__scan_stars__()

    def __scan_planets__(self) -> list[Planet]:
        """
        Returns a list of 9 visible planets [Sun...Pluto]
        """
        result = []
        for planet in range(10):
            raw_data = swe.calc_ut(self.jday, planet)
            result.append(Planet(
                lon=float(raw_data[0][0]),
                lat=float(raw_data[0][1])
            ))
        return result

    def __scan_stars__(self) -> list[Star]:
        """
        Returns a list of fixed stars
        """
        result = []
        for name in ALL_STARS:
            star = swe.fixstar_ut(
                name,
                self.jday
            )
            lon, lat = star[0][:2]
            result.append(Star(
                name=star[1],
                lon=lon,
                lat=lat
            ))
        return result


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':
    from datetime import datetime

    # Initiate a celestial sphere object
    test_sphere = Sphere(
        datetime(1976, 11, 29, 12, 15),
        time_zone=3,
        geo_lon=37 + 35/60,
        geo_lat=55 + 45/60,
    )
    stellar_objects = StellarObject(test_sphere)
    print(stellar_objects.planets[0])
    print(stellar_objects.stars[0])
