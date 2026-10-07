import numpy as np
import matplotlib.pyplot as plt
from Polarisability_Params import params
from scipy.constants import h, c, physical_constants


Eh = physical_constants['Hartree energy'][0]


def polarisability(np12_ns, np32_ns, dE12, dE32, A_ns, wavelength):

    omega = (h*c/wavelength)/Eh

    pol = 1/3*(
                ((((np12_ns)**2)*dE12)
               /
               ((dE12**2)-omega**2))

               +

               ((((np32_ns)**2)*dE32)
               /
               ((dE32**2)-omega**2))
               ) + A_ns

    return pol

 #species: n, np_1/2||D||ns, err np_1/2||D||ns,  np_3/2||D||ns, err np_3/2||D||ns, dE_(np_1/2), dE_(np_3/2), A_ns, dA_ns


wavelengths = np.arange(700e-9,1001e-9,1e-9)

species_list = [#'Li',
                #'Na',
                #'K',
                'Rb',
                'Cs',
                #'Fr'
                ]

for species in species_list:

    p = params[species]

    polaris = np.array(
        [
            polarisability(p[1], p[3], p[5], p[6], p[7], wavelength
            )
            for wavelength in wavelengths
        ]
    )

    polaris_low = np.array(
        [
            polarisability(p[1]+p[2], p[3]+p[4], p[5], p[6], p[7]+p[8], wavelength
            )
            for wavelength in wavelengths
            
        ]
    )

    polaris_high = np.array(
        [
            polarisability(p[1]-p[2], p[3]-p[4], p[5], p[6], p[7]-p[8], wavelength
            )
            for wavelength in wavelengths
            
        ]
    )

    err_polaris = np.abs(np.subtract(polaris_low,polaris))
    err_polaris_ = np.abs(np.subtract(polaris,polaris_high))

    plt.errorbar(wavelengths*1e9,
                 polaris,
                 (err_polaris,err_polaris_),
                 marker='x',
                 label=species)
    
plt.axhline(y=0, ls='--', linewidth = 0.5, color='black')

plt.xlabel('Tweezer wavelength (nm)')
plt.ylabel(r'Polarisability ($a_0^3$)')

plt.ylim([-5000,5000])
plt.legend()
plt.show()