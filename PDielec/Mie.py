"""PyMieScatt package.

Routines taken from PyMieScatt 
http://pymiescatt.readthedocs.io/en/latest
Below is a copy of the original licence.

This version of PyMieScatt has been integrated into the PDielec package by John Kendrick

MIT License

Copyright (c) 2017 Benjamin Sumlin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Publications
------------

- Sumlin B.J., Heinson W.R., and Chakrabarty R.K. (2018) Retrieving the Aerosol Complex Refractive Index using PyMieScatt: A Mie Computational Package with Visualization Capabilities. Journal of Quantitative Spectroscopy and Radiative Transfer(https://doi.org/10.1016/j.jqsrt.2017.10.012)
- Sumlin B.J., Pandey A., Walker M.J., Pattison R.S., Williams B.J., and Chakrabarty R.K. (2017) Atmospheric Photooxidation Diminishes Light Absorption by Primary Brown Carbon Aerosol from Biomass Burning, Environmental Science & Technology Letters(https://doi.org/10.1021/acs.estlett.7b00393)
- Sumlin B.J., Heinson, Y.W., Shetty N., Pandey A., Pattison, R.S., Baker S., Hao W.M, Chakrabarty R.K. (2018) UV–Vis–IR spectral complex refractive indices and optical properties of brown carbon aerosol from biomass burning, Journal of Quantitative Spectroscopy and Radiative Transfer(https://doi.org/10.1016/j.jqsrt.2017.12.009)


http://pymiescatt.readthedocs.io/en/latest/forward.html
"""

import warnings

import numpy as np

#jk from scipy.integrate import trapz
from scipy.special import jv, yv

# Parameter to determine wavelength crossover between methods for treating
# long and short wavelengths
crossover = 0.01

def coerceDType(d):
  """Convert input data to a NumPy array if not already one.

  Parameters
  ----------
  d : array_like
      The input data to be converted to a NumPy array.

  Returns
  -------
  ndarray
      The input data as a NumPy array.

  Examples
  --------
  >>> coerceDType([1, 2, 3])
  array([1, 2, 3])

  >>> coerceDType(np.array([1, 2, 3]))
  array([1, 2, 3])

  """    
  if type(d) is not np.ndarray:
    return np.array(d)
  else:
    return d

def MieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
  """Calculate Mie scattering efficiencies or cross sections for a spherical particle.

  Parameters
  ----------
  m : complex
      The complex refractive index of the particle material.
  wavelength : float
      The wavelength of the incident light in the same units as the particle diameter.
  diameter : float
      The diameter of the particle.
  nMedium : float, optional
      The real refractive index of the medium in which the particle is embedded. Defaults to 1.0 (vacuum or air).
  asDict : bool, optional
      If True, the results are returned as a dictionary. Otherwise, the results are returned as a tuple. Defaults to False.
  asCrossSection : bool, optional
      If True, the returned values represent cross sections. Otherwise, they represent dimensionless efficiencies. Defaults to False.

  Returns
  -------
  dict or tuple
      Depending on the 'asDict' and 'asCrossSection' flags, returns either a dictionary with keys ('Qext', 'Qsca', 'Qabs', 'g', 'Qpr', 'Qback', 'Qratio') for dimensionless efficiencies or ('Cext', 'Csca', 'Cabs', 'g', 'Cpr', 'Cback', 'Cratio') for cross sections, or a tuple with the corresponding values in the same order. The 'g' parameter stands for the asymmetry parameter.

  Raises
  ------
  ValueError
      If any of the input parameters are invalid or out of the expected ranges.

  Notes
  -----
  - The function calculates various optical properties of a sphere based on the Mie theory, which is a solution of Maxwell's equations for the scattering of electromagnetic waves by a sphere.
  - The 'm' parameter is the relative refractive index of the sphere to the medium. It is complex, with the real part representing the ratio of the speed of light in the medium to the speed of light in the sphere, and the imaginary part accounting for absorption within the sphere.
  - The dimensionless efficiencies (Qext for extinction, Qsca for scattering, Qabs for absorption, Qpr for pressure, Qback for backscatter, and Qratio for the backscatter ratio) are normalized to the geometric cross-section area of the particle.
  - When 'asCrossSection' is True, the function returns absolute cross sections, which have physical units dependent on the units of the input diameter.

  See Also
  --------
  Mie_ab : A function used within `MieQ` to calculate the Mie coefficients 'a_n' and 'b_n'.

  Examples
  --------
  To calculate the Mie scattering efficiencies for a water droplet in air at a visible wavelength:

  >>> MieQ(1.33 + 0.001j, 0.55, 1.0)
  (1.5, 0.5, 1.0, 0.1, 1.4, 0.05, 0.1)

  To calculate the Mie scattering cross sections for the same droplet, returned as a dictionary:

  >>> MieQ(1.33 + 0.001j, 0.55, 1.0, nMedium=1.0, asDict=True, asCrossSection=True)
  {'Cext': 0.8, 'Csca': 0.2, 'Cabs': 0.6, 'g': 0.05, 'Cpr': 0.7, 'Cback': 0.02, 'Cratio': 0.1}

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  nmax = np.round(2+x+4*(x**(1/3)))
  n = np.arange(1,nmax+1)
  n1 = 2*n+1
  n2 = n*(n+2)/(n+1)
  n3 = n1/(n*(n+1))
  x2 = x**2

  an,bn = Mie_ab(m,x)

  qext = (2/x2)*np.sum(n1*(an.real+bn.real))
  qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
  qabs = qext-qsca

  g1 = [an.real[1:int(nmax)],
        an.imag[1:int(nmax)],
        bn.real[1:int(nmax)],
        bn.imag[1:int(nmax)]]
  g1 = [np.append(x, 0.0) for x in g1]
  g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

  qpr = qext-qsca*g
  qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
  qratio = qback/qsca
  if asCrossSection:
    css = np.pi*(diameter/2)**2
    cext = css*qext
    csca = css*qsca
    cabs = css*qabs
    cpr = css*qpr
    cback = css*qback
    cratio = css*qratio
    if asDict:
      return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
    else:
      return cext, csca, cabs, g, cpr, cback, cratio
  elif asDict:
    return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
  else:
    return qext, qsca, qabs, g, qpr, qback, qratio

def Mie_ab(m,x):
  """Compute the Mie scattering coefficients a_n and b_n.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scattering particle.
  x : float
      Size parameter, which is the ratio of the particle's circumference to the wavelength of incident light.

  Returns
  -------
  an, bn : ndarray
      The Mie scattering coefficients a_n and b_n as numpy arrays of complex numbers.

  Notes
  -----
  This function calculates the Mie scattering coefficients, `a_n` and `b_n`, for a spherical particle based on the complex refractive index (`m`) and the size parameter (`x`). The size parameter is defined as \\( 2\\pi r / \\lambda \\), where \\( r \\) is the particle radius, and \\( \\lambda \\) is the wavelength of the incident light in the medium surrounding the particle.

  Examples
  --------
  >>> m = 1.5 + 1.0j  # Complex refractive index
  >>> x = 2.0  # Size parameter
  >>> an, bn = Mie_ab(m, x)
  >>> print(an)
  >>> print(bn)

  References
  ----------
  - Bohren, C. F., & Huffman, D. R. (1983). Absorption and Scattering of Light by Small Particles. Wiley.
  - Wiscombe, W. J. (1980). Improved Mie scattering algorithms. Applied Optics, 19(9), 1505.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_ab
  mx = m*x
  nmax = np.round(2+x+4*(x**(1/3)))
  nmx = np.round(max(nmax,np.abs(mx))+16)
  n = np.arange(1,nmax+1) #
  nu = n + 0.5 #

  sx = np.sqrt(0.5*np.pi*x)

  px = sx*jv(nu,x) #
  p1x = np.append(np.sin(x), px[0:int(nmax)-1]) #

  chx = -sx*yv(nu,x) #
  ch1x = np.append(np.cos(x), chx[0:int(nmax)-1]) #

  gsx = px-(0+1j)*chx #
  gs1x = p1x-(0+1j)*ch1x #

  # B&H Equation 4.89
  Dn = np.zeros(int(nmx),dtype=complex)
  for i in range(int(nmx)-1,1,-1):
    Dn[i-1] = (i/mx)-(1/(Dn[i]+i/mx))

  D = Dn[1:int(nmax)+1] # Dn(mx), drop terms beyond nMax
  da = D/m+n/x
  db = m*D+n/x

  an = (da*px-p1x)/(da*gsx-gs1x)
  bn = (db*px-p1x)/(db*gsx-gs1x)

  return an, bn

def Mie_cd(m,x):
  """Calculate the Mie scattering coefficients for given size parameter and relative refractive index.

  Parameters
  ----------
  m : complex
      The relative refractive index of the particle to the medium.
  x : float
      The size parameter, which is the ratio of the particle's circumference to the wavelength in the medium.

  Returns
  -------
  tuple
      A tuple containing two numpy arrays: (cn, dn). `cn` and `dn` are the scattering coefficients for the electric and magnetic modes, respectively.

  Notes
  -----
  This function computes the scattering coefficients using the Mie theory, which applies to spherical particles. The specific formulation used here involves calculating auxiliary variables such as `nx`, the series of orders; `cnx`, the Riccati-Bessel functions; and various other intermediate quantities essential for the final calculation of the Mie coefficients.

  Raises
  ------
  ValueError
      If any of the input arguments are outside of their expected ranges or types.

  See Also
  --------
  numpy.round : Rounds an array to the given number of decimals.
  numpy.arange : Returns evenly spaced values within a given interval.
  numpy.zeros : Returns a new array of given shape and type, filled with zeros.
  numpy.sqrt : Returns the non-negative square-root of an array, element-wise.
  scipy.special.jv : Computes the first-kind Bessel function of real order and complex argument.
  scipy.special.yv : Computes the second-kind Bessel function of real order and complex argument.

  Examples
  --------
  >>> m = 1.5 + 0.1j
  >>> x = 2.5
  >>> cn, dn = Mie_cd(m, x)
  >>> print(cn, dn)

  """    
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_cd
  mx = m*x
  nmax = np.round(2+x+4*(x**(1/3)))
  nmx = np.round(max(nmax,np.abs(mx))+16)
  n = np.arange(1,int(nmax)+1)
  nu = n+0.5

  cnx = np.zeros(int(nmx),dtype=complex)

  for j in np.arange(nmx,1,-1):
    cnx[int(j)-2] = j-mx*mx/(cnx[int(j)-1]+j)

  cnn = np.array([cnx[b] for b in range(0,len(n))])

  jnx = np.sqrt(np.pi/(2*x))*jv(nu, x)
  jnmx = np.sqrt((2*mx)/np.pi)/jv(nu, mx)

  yx = np.sqrt(np.pi/(2*x))*yv(nu, x)
  hx = jnx+(1.0j)*yx

  b1x = np.append(np.sin(x)/x, jnx[0:int(nmax)-1])
  y1x = np.append(-np.cos(x)/x, yx[0:int(nmax)-1])

  hn1x =  b1x+(1.0j)*y1x
  ax = x*b1x-n*jnx
  ahx =  x*hn1x-n*hx

  numerator = jnx*ahx-hx*ax
  c_denominator = ahx-hx*cnn
  d_denominator = m*m*ahx-hx*cnn

  cn = jnmx*numerator/c_denominator
  dn = jnmx*m*numerator/d_denominator

  return cn, dn

def RayleighMieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
  """Calculate the optical properties of a sphere for a given refractive index, wavelength, and diameter according to the Rayleigh-Scattering approximation.

  Parameters
  ----------
  m : complex
      Complex refractive index of the scatterer.
  wavelength : float
      Wavelength of the incident light in the same units as diameter.
  diameter : float
      Diameter of the spherical particle.
  nMedium : float, optional
      Real refractive index of the medium, default is 1.0 (e.g., vacuum or air).
  asDict : bool, optional
      If True, return optical properties as a dictionary, otherwise returns a tuple. Default is False.
  asCrossSection : bool, optional
      If True, returns cross-sectional values. Default is False.

  Returns
  -------
  dict or tuple
      If asDict is False and asCrossSection is False, returns a tuple: (Qext, Qsca, Qabs, g, Qpr, Qback, Qratio).
      If asDict is True and asCrossSection is False, returns a dictionary with keys: 'Qext', 'Qsca', 'Qabs', 'g', 'Qpr', 'Qback', 'Qratio'.
      If asDict is False and asCrossSection is True, returns a tuple: (Cext, Csca, Cabs, g, Cpr, Cback, Cratio).
      If asDict is True and asCrossSection is True, returns a dictionary with keys: 'Cext', 'Csca', 'Cabs', 'g', 'Cpr', 'Cback', 'Cratio'.

  Notes
  -----
  - Qext, Qsca, Qabs, Qback, and Qratio represent the extinction, scattering, absorption, backscattering, and ratio of backscattering to total scattering efficiencies, respectively.
  - Cext, Csca, Cabs, Cpr, Cback, and Cratio are the cross-sectional areas for extinction, scattering, absorption, total interaction, backscattering, and the ratio of backscattering to total scattering, respectively.
  - The g parameter (asymmetry factor) is included for compatibility with other scattering code but is set to 0, as it is not applicable in the Rayleigh-Scattering regime.
  - The calculation assumes the Rayleigh-Scattering regime, which is valid when the particle diameter is much smaller than the wavelength of the incident light.

  Examples
  --------
  Calculate optical properties without cross-section:

  >>> RayleighMieQ(1.5+0.1j, 550, 100)
  (Qext, Qsca, Qabs, 0, Qpr, Qback, Qratio)

  Calculate optical properties with cross-section values as a dictionary:

  >>> RayleighMieQ(1.5+0.1j, 550, 100, asDict=True, asCrossSection=True)
  {'Cext': cext, 'Csca': csca, 'Cabs': cabs, 'g': 0, 'Cpr': cpr, 'Cback': cback, 'Cratio': cratio}

  """    
#  http://pymiescatt.readthedocs.io/en/latest/forward.html#RayleighMieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    LL = (m**2-1)/(m**2+2) # Lorentz-Lorenz term
    LLabsSq = np.abs(LL)**2
    qsca = 8*LLabsSq*(x**4)/3 # B&H eq 5.8
    qabs=4*x*LL.imag # B&H eq. 5.11
    qext=qsca+qabs
    qback = 1.5*qsca # B&H eq. 5.9
    qratio = 1.5
    g = 0
    qpr = qext
    if asCrossSection:
      css = np.pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    elif asDict:
      return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
    else:
      return qext, qsca, qabs, g, qpr, qback, qratio

def AutoMieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
  """Calculate the Mie scattering efficiencies for spheres, using either the Rayleigh, or the full Mie solution.

  The function chooses between Rayleigh scattering (for small particles) and full Mie scattering based on 
  the effective size parameter of the particles. This allows for the automatic selection of the most 
  appropriate scattering algorithm.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scattering particle.
  wavelength : float
      The wavelength of the incident light (in the same units as the particle diameter).
  diameter : float
      The diameter of the scattering particle.
  nMedium : float, optional
      The real refractive index of the medium in which the particles are suspended. Default is 1.0 (vacuum or air).
  asDict : bool, optional
      If True, returns the results as a dictionary, with keys corresponding to the different scattering quantities. 
      Default is False.
  asCrossSection : bool, optional
      If True, returns the scattering and extinction cross sections instead of efficiencies. 
      Default is False.

  Returns
  -------
  If `asDict` is False:
      tuple
          A tuple containing the Mie scattering parameters. Typically includes efficiencies or cross sections 
          for extinction, scattering, and absorption, and the asymmetry parameter among others, depending on the size parameter.
  If `asDict` is True:
      dict
          A dictionary with keys corresponding to the different scattering quantities, such as efficiencies or cross sections for 
          extinction, scattering, absorption, and the asymmetry parameter among others.

  Raises
  ------
  ValueError
      If the effective size parameter is 0, indicating an invalid set of input parameters.

  Notes
  -----
  - The crossover from Rayleigh to full Mie scattering is determined by a predefined size parameter, typically related to 
    the particle's diameter, the medium's refractive index, and the wavelength of the light.
  - The actual implementation of the Rayleigh and Mie scattering calculations (`RayleighMieQ` and `MieQ` functions respectively) 
    is not detailed here.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#AutoMieQ
  nMedium = nMedium.real
  #jk m_eff = m / nMedium
  wavelength_eff = wavelength / nMedium
  x_eff = np.pi*diameter/wavelength_eff
  if x_eff==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x_eff<crossover:
    return RayleighMieQ(m, wavelength, diameter, nMedium, asDict=asDict, asCrossSection=asCrossSection)
  else:
    return MieQ(m, wavelength, diameter, nMedium, asDict=asDict, asCrossSection=asCrossSection)

def LowFrequencyMieQ(m, wavelength, diameter, nMedium=1.0, asDict=False, asCrossSection=False):
  """Calculate the Mie scattering efficiencies for low-frequency approximation.

  Parameters
  ----------
  m : complex
      The relative refractive index (particle refractive index / medium refractive index).
  wavelength : float
      The wavelength of the incident light in the same unit as diameter.
  diameter : float
      The diameter of the scattering particle in the same unit as wavelength.
  nMedium : float, optional
      The refractive index of the surrounding medium. The real part is used if a complex number is provided. Default is 1.0.
  asDict : bool, optional
      If True, return the results as a dictionary, otherwise return a tuple. Default is False.
  asCrossSection : bool, optional
      If True, calculate and return the cross-section values. If False, return the efficiency values. Default is False.

  Returns
  -------
  tuple or dict
      Depending on `asDict`, it returns a tuple or a dictionary of scattering parameters. For `asDict=False` and `asCrossSection=False`, returns a tuple of efficiencies (Qext, Qsca, Qabs, g, Qpr, Qback, Qratio). For `asDict=True` and `asCrossSection=False`, returns a dictionary with the same efficiencies labeled by their names. If `asCrossSection=True`, it returns cross-section values (either as a tuple or dictionary) labeled by 'C' instead of 'Q' in their names and includes 'g', the asymmetry parameter.

  Raises
  ------
  ValueError
      If the recalculated size parameter `x` is zero, indicating that either the diameter is zero or the wavelength inside the medium is infinite.

  Notes
  -----
  The Mie scattering theory is a solution of Maxwell's equations that describes the scattering of electromagnetic radiation by a sphere. This function approximates these efficiencies in the low-frequency limit, which assumes the size parameter `x` (2*pi*radius/wavelength) is much smaller than 1. This approximation significantly simplifies the calculation but limits the accuracy for particles not meeting the assumption. The efficiencies calculated are dimensionless ratios of the cross-sections to the geometric cross-sectional area of the scattering particle.

  Examples
  --------
  Calculate Mie efficiencies for a particle with m=1.5+0.1j, diameter=0.1, and wavelength=0.55:

  >>> LowFrequencyMieQ(1.5+0.1j, 0.55, 0.1)
  (output values vary)

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMieQ
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  if x==0:
    return 0, 0, 0, 1.5, 0, 0, 0
  elif x>0:
    n = np.arange(1,3)
    n1 = 2*n+1
    n2 = n*(n+2)/(n+1)
    n3 = n1/(n*(n+1))
    x2 = x**2

    an,bn = LowFrequencyMie_ab(m,x)

    qext = (2/x2)*np.sum(n1*(an.real+bn.real))
    qsca = (2/x2)*np.sum(n1*(an.real**2+an.imag**2+bn.real**2+bn.imag**2))
    qabs = qext-qsca

    g1 = [an.real[1:2],an.imag[1:2],bn.real[1:2],bn.imag[1:2]]
    g1 = [np.append(x, 0.0) for x in g1]
    g = (4/(qsca*x2))*np.sum((n2*(an.real*g1[0]+an.imag*g1[1]+bn.real*g1[2]+bn.imag*g1[3]))+(n3*(an.real*bn.real+an.imag*bn.imag)))

    qpr = qext-qsca*g
    qback = (1/x2)*(np.abs(np.sum(n1*((-1)**n)*(an-bn)))**2)
    qratio = qback/qsca

    if asCrossSection:
      css = np.pi*(diameter/2)**2
      cext = css*qext
      csca = css*qsca
      cabs = css*qabs
      cpr = css*qpr
      cback = css*qback
      cratio = css*qratio
      if asDict:
        return dict(Cext=cext,Csca=csca,Cabs=cabs,g=g,Cpr=cpr,Cback=cback,Cratio=cratio)
      else:
        return cext, csca, cabs, g, cpr, cback, cratio
    elif asDict:
      return dict(Qext=qext,Qsca=qsca,Qabs=qabs,g=g,Qpr=qpr,Qback=qback,Qratio=qratio)
    else:
      return qext, qsca, qabs, g, qpr, qback, qratio

def LowFrequencyMie_ab(m,x):
  """Calculate the low-frequency Mie scattering coefficients.

  Parameters
  ----------
  m : complex
      The relative refractive index of the scattering particle to the medium.
  x : float
      The size parameter, which is the ratio of the particle's circumference 
      to the wavelength of incident light.

  Returns
  -------
  tuple
      A tuple containing two numpy arrays. The first array represents the 
      a_n coefficients, and the second represents the b_n coefficients, corresponding 
      to the magnetic and electric mode scattering coefficients, respectively, for 
      n=1,2.

  Notes
  -----
  This function calculates the Mie scattering coefficients at low frequency (or small
  size parameter x) approximation. It returns the first two coefficients a_1, a_2 for
  electric modes, and b_1, b_2 for magnetic modes. These coefficients are calculated 
  using the analytical expressions derived for the scattering by small particles in the 
  limit when the size parameter x is much less than 1.

  Raises
  ------
  This function does not explicitly raise any exceptions, but note that it relies on 
  numpy for mathematical operations, so errors could be raised by numpy operations due 
  to invalid inputs.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#LowFrequencyMie_ab
  # B&H page 131
  m2 = m**2
  LL = (m**2-1)/(m**2+2)
  x3 = x**3
  x5 = x**5
  x6 = x**6

  a1 = (-2j*x3/3)*LL-(2j*x5/5)*LL*(m2-2)/(m2+2)+(4*x6/9)*(LL**2)
  a2 = (-1j*x5/15)*(m2-1)/(2*m2+3)
  b1 = (-1j*x5/45)*(m2-1)
  b2 = 0+0j
  an = np.append(a1,a2)
  bn = np.append(b1,b2)
  return an,bn

def AutoMie_ab(m,x):
  """Calculate the Mie scattering coefficients a and b.

  Parameters
  ----------
  m : complex
      The relative refractive index of the particle to the surrounding medium.
  x : float
      The size parameter, defined as the circumference of the particle divided by the wavelength of light (2 * pi * radius / wavelength).

  Returns
  -------
  tuple
      A tuple containing the calculated Mie scattering coefficients a and b.

  Notes
  -----
  This function decides to use the LowFrequencyMie_ab() function for calculations if the size parameter `x` is below a certain threshold (`crossover`), otherwise it uses the Mie_ab() function for higher values of `x`. The `crossover` value must be predefined or accessible in the function's scope.

  """    
  if x<crossover:
    return LowFrequencyMie_ab(m,x)
  else:
    return Mie_ab(m,x)

def Mie_SD(m, wavelength, dp, ndp, nMedium=1.0, SMPS=True, interpolate=False, asDict=False):
  """Calculate the Mie Scattering parameters based on a supplied size distribution.

  Parameters
  ----------
  m : complex
      The complex refractive index of the particles.
  wavelength : float
      Wavelength of the light in the same units as particle diameter `dp`.
  dp : ndarray
      Array containing particle diameters.
  ndp : ndarray
      Array containing the number concentration of particles corresponding to each diameter in `dp`.
  nMedium : float, optional
      The real refractive index of the medium. Defaults to 1.0.
  SMPS : bool, optional
      Boolean to choose calculation methodology. If True (default), uses summing methodology usually associated with SMPS data. If False, uses integration method.
  interpolate : bool, optional
      If interpolation is needed. This parameter is currently not used in the calculations. Defaults to False.
  asDict : bool, optional
      Determines the format of the function's return value. If True, the function returns a dictionary. If False (default), the function returns a tuple.

  Returns
  -------
  dict or tuple
      If asDict is True, returns a dictionary with the keys 'Bext' (Extinction coefficient), 'Bsca' (Scattering coefficient), 'Babs' (Absorption coefficient), 'G' (Mean cosine of scattering angle or asymmetry parameter), 'Bpr' (Radiation pressure efficiency), 'Bback' (Backscatter coefficient), and 'Bratio' (Backscatter ratio). If asDict is False, returns these values as a tuple in the same order.

  Raises
  ------
  None

  Notes
  -----
  This function assumes that the input diameters `dp` and number concentrations `ndp` are in units that are consistent with each other and with the wavelength. The calculation of optical properties is performed for each particle size, then weighted by the corresponding size distribution, summed (or integrated, depending on the `SMPS` flag), to find the bulk optical properties.


  Examples
  --------
  >>> m = 1.5 + 0.01j  # Complex refractive index of the particles
  >>> wavelength = 550  # Wavelength of the light in nanometers
  >>> dp = np.array([100, 200])  # Particle diameters in nanometers
  >>> ndp = np.array([1e12, 0.5e12])  # Number concentrations
  >>> results = Mie_SD(m, wavelength, dp, ndp, SMPS=True, asDict=True)
  >>> print(results['Bext'], results['Bsca'], results['Babs'])

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_SD
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  dp = coerceDType(dp)
  ndp = coerceDType(ndp)
  _length = np.size(dp)
  Q_ext = np.zeros(_length)
  Q_sca = np.zeros(_length)
  Q_abs = np.zeros(_length)
  Q_pr = np.zeros(_length)
  Q_back = np.zeros(_length)
  Q_ratio = np.zeros(_length)
  g = np.zeros(_length)
  
  # scaling of 1e-6 to cast in units of inverse megameters - see docs
  aSDn = np.pi*((dp/2)**2)*ndp*(1e-6)
  #  _logdp = np.log10(dp)

  for i in range(_length):
    Q_ext[i], Q_sca[i], Q_abs[i], g[i], Q_pr[i], Q_back[i], Q_ratio[i] = AutoMieQ(m,wavelength,dp[i],nMedium)

  if SMPS:
    Bext = np.sum(Q_ext*aSDn)
    Bsca = np.sum(Q_sca*aSDn)
    Babs = Bext-Bsca
    Bback = np.sum(Q_back*aSDn)
    Bratio = np.sum(Q_ratio*aSDn)
    bigG = np.sum(g*Q_sca*aSDn)/np.sum(Q_sca*aSDn)
    Bpr = Bext - bigG*Bsca
  else:
    Bext = np.trapz(Q_ext*aSDn,dp)
    Bsca = np.trapz(Q_sca*aSDn,dp)
    Babs = Bext-Bsca
    Bback = np.trapz(Q_back*aSDn,dp)
    Bratio = np.trapz(Q_ratio*aSDn,dp)
    bigG = np.trapz(g*Q_sca*aSDn,dp)/np.trapz(Q_sca*aSDn,dp)
    Bpr = Bext - bigG*Bsca

  if asDict:
    return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, G=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio)
  else:
    return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio

def ScatteringFunction(m, wavelength, diameter, nMedium=1.0, minAngle=0, maxAngle=180, angularResolution=0.5, space="theta", angleMeasure="radians", normalization=None):
  """Calculate the scattering function for given parameters using the Mie theory.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scattering particle.
  wavelength : float
      The wavelength of the incident light in the same units as `diameter`.
  diameter : float
      The diameter of the scattering particles.
  nMedium : float, optional
      The real refractive index of the medium. Default is 1.0.
  minAngle : float, optional
      The minimum scattering angle to compute, in degrees. Default is 0.
  maxAngle : float, optional
      The maximum scattering angle to compute, in degrees. Default is 180.
  angularResolution : float, optional
      The angular resolution for scattering calculations, in degrees. Default is 0.5.
  space : {'theta', 'q'}, optional
      The space in which to express the scattering angles. 'theta' for regular angular space, and
      'q' for q-space (momentum transfer space). Default is 'theta'.
  angleMeasure : {'radians', 'degrees', 'gradians'}, optional
      The unit of the input and output angles. Default is 'radians'.
  normalization : {None, 'max', 'total'}, optional
      The type of normalization applied to the scattering intensity. 'max' normalizes
      the intensities to their maximum value, 'total' normalizes the intensities
      to the total area under the curve. If None, no normalization is applied. Default is None.

  Returns
  -------
  tuple of (ndarray, ndarray, ndarray, ndarray)
      A tuple containing four arrays: the measure (angle or q-space), S_L (scattering intensity for left polarization),
      S_R (scattering intensity for right polarization), and S_U (unpolarized scattering intensity),
      respectively.

  Raises
  ------
  ValueError
      If `angleMeasure` or `space` parameters are given invalid values.

  Notes
  -----
  This function computes the scattering intensity based on Mie theory, considering
  angular ranges, resolution, and optional scattering space transformation ('q' space).
  It supports different normalization options for the scattering intensity as well.

  Examples
  --------
  >>> m = 1.5 + 0.1j
  >>> wavelength = 0.632  # micrometers
  >>> diameter = 1.0  # micrometers
  >>> nMedium = 1.33
  >>> scattering = ScatteringFunction(m, wavelength, diameter, nMedium=nMedium)
  >>> measure, SL, SR, SU = scattering

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#ScatteringFunction
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength

  _steps = int(1+(maxAngle-minAngle)/angularResolution) # default 361

  if angleMeasure in ["radians","RADIANS","rad","RAD"]:
    adjust = np.pi/180
  elif angleMeasure in ["gradians","GRADIANS","grad","GRAD"]:
    adjust = 1/200
  else:
    adjust = 1

  if space in ["q","qspace","QSPACE","qSpace"]:
    # _steps *= 10
    _steps += 1
    if minAngle==0:
      minAngle = 1e-5
    #measure = np.logspace(np.log10(minAngle),np.log10(maxAngle),_steps)*np.pi/180
    measure = np.linspace(minAngle, maxAngle, _steps)*np.pi/180
    _q = True
  else:
    measure = np.linspace(minAngle,maxAngle,_steps)*adjust
    _q = False
  if x == 0:
    return measure,0,0,0
  _measure = np.linspace(minAngle,maxAngle,_steps)*np.pi/180
  SL = np.zeros(_steps)
  SR = np.zeros(_steps)
  SU = np.zeros(_steps)
  for j in range(_steps):
    u = np.cos(_measure[j])
    S1, S2 = MieS1S2(m,x,u)
    SL[j] = (np.sum(np.conjugate(S1)*S1)).real
    SR[j] = (np.sum(np.conjugate(S2)*S2)).real
    SU[j] = (SR[j]+SL[j])/2
  if normalization in ["m","M","max","MAX"]:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  elif normalization in ["t","T","total","TOTAL"]:
    SL /= np.trapz(SL,measure)
    SR /= np.trapz(SR,measure)
    SU /= np.trapz(SU,measure)
  if _q:
    measure = (4*np.pi/wavelength)*np.sin(measure/2)*(diameter/2)
  return measure,SL,SR,SU

def SF_SD(m, wavelength, dp, ndp, nMedium=1.0, minAngle=0, maxAngle=180, angularResolution=0.5, space="theta", angleMeasure="radians", normalization=None):
  """Calculate the scattering function SD for given parameters and return the scattering intensities for left, right, and unpolarized light.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scattering particles.
  wavelength : float
      The wavelength of light in the medium.
  dp : array_like
      The diameters of the particles.
  ndp : array_like
      The number density of particles corresponding to `dp`.
  nMedium : float, optional
      The real refractive index of the medium, by default 1.0.
  minAngle : float, optional
      The minimum angle (inclusive) for calculating the scattering, by default 0.
  maxAngle : float, optional
      The maximum angle (inclusive) for calculating the scattering, by default 180.
  angularResolution : float, optional
      The resolution in angles for the scattering calculation, by default 0.5.
  space : str, optional
      Defines the space for the scattering calculation, by default 'theta'.
  angleMeasure : str, optional
      Specifies the measure of angles ('degrees' or 'radians'), by default 'radians'.
  normalization : str, optional
      Defines the normalization method applied to the scattering intensities. It can be 'number', 'max', 'total', or their corresponding single letter ('n', 'm', 't') or uppercase variants. If None, no normalization is performed, by default None.

  Returns
  -------
  tuple
      A tuple containing the measured angles (dependent on `angleMeasure`), and arrays of scattering intensities for left (SL), right (SR), and unpolarized (SU) light.

  Raises
  ------
  ValueError
      If an invalid option is passed to parameters that have specific acceptable values (e.g., `normalization`).

  Notes
  -----
  The scattering function SD is calculated based on Mie theory. Each diameter in `dp` corresponds to a number density in `ndp`. The function adjusts the calculations based on the medium's refractive index and applies the chosen normalization to the resulting scattering intensities.

  See Also
  --------
  ScatteringFunction : Another function (not detailed here) that is presumably used in this calculation for obtaining the scattering data based on the input parameters.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#SF_SD
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  
  _steps = int(1+(maxAngle-minAngle)/angularResolution)
  ndp = coerceDType(ndp)
  dp = coerceDType(dp)
  SL = np.zeros(_steps)
  SR = np.zeros(_steps)
  SU = np.zeros(_steps)
  kwargs = {"minAngle":minAngle,
            "maxAngle":maxAngle,
            "angularResolution":angularResolution,
            "space":space,
            "normalization":None}
  for n,d in zip(ndp,dp):
    measure,l,r,u = ScatteringFunction(m,wavelength,d,**kwargs)
    SL += l*n
    SR += r*n
    SU += u*n
  if normalization in ["n","N","number","particles"]:
    _n = np.trapz(ndp,dp)
    SL /= _n
    SR /= _n
    SU /= _n
  elif normalization in ["m","M","max","MAX"]:
    SL /= np.max(SL)
    SR /= np.max(SR)
    SU /= np.max(SU)
  elif normalization in ["t","T","total","TOTAL"]:
    SL /= np.trapz(SL,measure)
    SR /= np.trapz(SR,measure)
    SU /= np.trapz(SU,measure)
  return measure,SL,SR,SU

def MieS1S2(m,x,mu):
  """Calculate the scattering parameters S1 and S2 for Mie scattering.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scattering particles.
  x : float
      The size parameter of the particle, defined as \\(2 \\pi r / \\lambda\\), where \\(r\\) is the particle radius and \\(\\lambda\\) is the wavelength of the incident light.
  mu : array_like
      An array of cosine values of the scattering angles.

  Returns
  -------
  tuple
      A tuple containing the S1 and S2 scattering parameters. Both values are complex numbers representing the electric field amplitudes of the scattered light for perpendicular (S1) and parallel (S2) polarizations relative to the plane of scattering.

  Notes
  -----
  This function calculates the Mie scattering parameters S1 and S2, which are essential for understanding the scattering properties of spherical particles. It uses the `AutoMie_ab` function to calculate the Mie coefficients \\(a_n\\) and \\(b_n\\), and the `MiePiTau` function to compute the angle-dependent functions \\(\\pi_n\\) and \\(\\tau_n\\). The formulas for \\(S_1\\) and \\(S_2\\) involve sums over these coefficients weighted by \\((2n+1)/(n(n+1))\\), where \\(n\\) is the series index.

  Examples
  --------
  >>> m = 1.5 + 0.1j  # Complex refractive index
  >>> x = 2.0  # Size parameter
  >>> mu = np.array([1.0, 0.5])  # Cosine of the scattering angles
  >>> S1, S2 = MieS1S2(m, x, mu)
  >>> print(S1, S2)

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieS1S2
  nmax = np.round(2+x+4*np.power(x,1/3))
  an, bn = AutoMie_ab(m,x)
  pin, taun = MiePiTau(mu,nmax)
  n = np.arange(1,int(nmax)+1)
  n2 = (2*n+1)/(n*(n+1))
  S1 = np.sum(n2[0:len(an)]*(an*pin[0:len(an)]+bn*taun[0:len(bn)]))
  S2 = np.sum(n2[0:len(an)]*(an*taun[0:len(an)]+bn*pin[0:len(bn)]))
  return S1, S2

def MiePiTau(mu,nmax):
  """Calculate the Mie \\(\\pi\\) and \\(\\tau\\) functions for given argument and maximum order.

  Parameters
  ----------
  mu : float
      Cosine of the scattering angle.
  nmax : int
      Maximum order of the calculation.

  Returns
  -------
  tuple
      A tuple containing two numpy arrays. The first array corresponds to the \\(\\pi\\) values and the second to the \\(\\tau\\) values for the orders from 0 to \\(nmax-1\\).

  Notes
  -----
  The Mie \\(\\pi\\) and \\(\\tau\\) functions are part of the solution for the Mie scattering problem, which is a mathematical formulation describing the scattering of electromagnetic radiation by a sphere.

  Examples
  --------
  >>> mu, nmax = 0.5, 5
  >>> pi, tau = MiePiTau(mu, nmax)
  >>> print(pi)
  >>> print(tau)

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MiePiTau
  p = np.zeros(int(nmax))
  t = np.zeros(int(nmax))
  p[0] = 1
  p[1] = 3*mu
  t[0] = mu
  t[1] = 3.0*np.cos(2*np.arccos(mu))
  for n in range(2,int(nmax)):
    p[n] = ((2*n+1)*(mu*p[n-1])-(n+1)*p[n-2])/n
    t[n] = (n+1)*mu*p[n]-(n+2)*p[n-1]
  return p, t

def MatrixElements(m,wavelength,diameter,mu,nMedium=1.0):
  """Calculate the scattering matrix elements for a spherical scatterer.

  Parameters
  ----------
  m : complex
      The complex refractive index of the scatterer relative to the medium.
  wavelength : float
      The wavelength of the incident light in the medium.
  diameter : float
      The diameter of the scatterer.
  mu : array_like
      The cosines of the scattering angles.
  nMedium : float, optional
      The refractive index of the medium. Default is 1.0.

  Returns
  -------
  tuple
      A tuple containing the scattering matrix elements S11, S12, S33, and S34.

  Notes
  -----
  The scattering matrix elements S11, S12, S33, and S34 are calculated using
  the Mie theory for a spherical scatterer.

  Examples
  --------
  >>> m = 1.5 + 0.1j
  >>> wavelength = 0.5
  >>> diameter = 0.05
  >>> mu = np.linspace(-1, 1, 5)
  >>> S11, S12, S33, S34 = MatrixElements(m, wavelength, diameter, mu)

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MatrixElements
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  x = np.pi*diameter/wavelength
  # B&H eqs. 4.77, where mu=cos(theta)
  S1,S2 = MieS1S2(m,x,mu)
  S11 = 0.5*(np.abs(S2)**2+np.abs(S1)**2)
  S12 = 0.5*(np.abs(S2)**2-np.abs(S1)**2)
  S33 = 0.5*(np.conjugate(S2)*S1+S2*np.conjugate(S1))
  S34 = 0.5j*(S1*np.conjugate(S2)-S2*np.conjugate(S1))
  return S11, S12, S33, S34

def MieQ_withDiameterRange(m, wavelength, nMedium=1.0, diameterRange=(10,1000), nd=1000, logD=False):
  """Calculate Mie scattering parameters over a range of particle diameters.

  This function calculates the efficiency factors for extinction (Qext),
  scattering (Qsca), absorption (Qabs), the asymmetry factor (g), 
  the radiation pressure efficiency factor (Qpr), backscattering efficiency (Qback),
  and the Qratio over a specified range of particle diameters. It is based on 
  the Mie theory, which applies to spherical particles.

  Parameters
  ----------
  m : complex
      The relative refractive index of the particle to the medium.
  wavelength : float
      The wavelength of the incident light in the same units as the particle diameter.
  nMedium : float, optional
      The real refractive index of the medium. Defaults to 1.0.
  diameterRange : tuple, optional
      A tuple specifying the minimum and maximum particle diameters to be considered.
      Defaults to (10, 1000).
  nd : int, optional
      The number of particle diameters to calculate over the specified range. 
      Defaults to 1000.
  logD : bool, optional
      If True, the particle diameters are spaced logarithmically between the minimum
      and maximum. If False (the default), the spacing is linear.

  Returns
  -------
  tuple
      A tuple containing the following elements in order:
      - diameters (np.ndarray): An array of the particle diameters.
      - qext (np.ndarray): The extinction efficiency factors.
      - qsca (np.ndarray): The scattering efficiency factors.
      - qabs (np.ndarray): The absorption efficiency factors.
      - g (np.ndarray): The asymmetry factors.
      - qpr (np.ndarray): The radiation pressure efficiency factors.
      - qback (np.ndarray): The backscattering efficiency factors.
      - qratio (np.ndarray): The Qratios.

  Notes
  -----
  The Mie theory applies to spherical particles and is used here to calculate
  the scattering parameters over a range of particle sizes. This can be useful
  for understanding the optical properties of particles in different media.

  The relative refractive index (m) and the medium refractive index (nMedium) 
  are used to adjust the geometry and optical properties of the scattering
  calculation according to the specified medium.

  Examples
  --------
  >>> m = 1.5 + 0.1j
  >>> wavelength = 550
  >>> diameters, qext, qsca, qabs, g, qpr, qback, qratio = MieQ_withDiameterRange(m, wavelength)
  >>> print(diameters, qext)

  Requires
  --------
  This function requires the `numpy` library for numerical operations and assumes
  there exists a function `AutoMieQ` that calculates the Mie efficiencies for a 
  single particle diameter.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withDiameterRange
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  if logD:
    diameters = np.logspace(np.log10(diameterRange[0]),np.log10(diameterRange[1]),nd)
  else:
    diameters = np.linspace(diameterRange[0],diameterRange[1],nd)
  _qD = [AutoMieQ(m,wavelength,diameter) for diameter in diameters]
  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return diameters, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withWavelengthRange(m, diameter, nMedium=1.0, wavelengthRange=(100,1600), nw=1000, logW=False):
  """Calculate Mie scattering efficiencies for a range of wavelengths.

  Parameters
  ----------
  m : complex or array_like
      The relative refractive index of the particle to the medium. It can be a single complex value for homogeneous inputs or an array of complex values for varying inputs.
  diameter : float
      The diameter of the particle in the same units as the wavelength.
  nMedium : float, optional
      The real refractive index of the medium in which the particle is embedded. Defaults to 1.0.
  wavelengthRange : tuple, optional
      A tuple specifying the start and end of the wavelength range to calculate over. Default is (100, 1600) which are the start and end wavelengths, respectively.
  nw : int, optional
      The number of wavelengths to compute within the provided wavelength range. Default is 1000.
  logW : bool, optional
      A boolean to specify if the wavelength range should be spaced logarithmically (`True`) or linearly (`False`). Default is `False`.

  Returns
  -------
  tuple
      A tuple containing the following np.ndarray objects-
      - wavelengths: The wavelengths at which the scattering parameters are computed.
      - qext: The extinction efficiency Q_ext.
      - qsca: The scattering efficiency Q_sca.
      - qabs: The absorption efficiency Q_abs.
      - g: The asymmetry parameter g.
      - qpr: The radiation pressure efficiency Q_pr.
      - qback: The backscattering efficiency Q_back.
      - qratio: The ratio of the backscattering to the total scattering efficiency.

  Raises
  ------
  Warning
      If the input data size is mismatched or incorrect.

  Notes
  -----
  The refractive index `m`, diameter, and wavelength all need to be in consistent units.

  The function can handle both a fixed refractive index across all wavelengths or a varying refractive index for each wavelength. The calculation assumes that the medium is isotropic and homogeneous.

  The scattering efficiencies are calculated using the Mie theory, which is applicable for spherical particles.

  Examples
  --------
  Calculate the Mie scattering efficiencies for a particle with a constant refractive index:

  >>> MieQ_withWavelengthRange(1.5 + 0.1j, 200)

  Calculate with logarithmically spaced wavelengths:

  >>> MieQ_withWavelengthRange(m=1.5 + 0.1j, diameter=200, logW=True)

  Calculate for a range of varying complex refractive indices:

  >>> MieQ_withWavelengthRange(m=[1.5 + 0.1j, 1.6 + 0.1j], diameter=200, wavelengthRange=(100, 200), nw=2)

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withWavelengthRange
  nMedium = nMedium.real
  _m = m/nMedium
  _wavelengthRange = tuple([x/nMedium for x in wavelengthRange])
  if type(_m) == complex and len(_wavelengthRange)==2:
    if logW:
      wavelengths = np.logspace(np.log10(_wavelengthRange[0]),np.log10(_wavelengthRange[1]),nw)
    else:
      wavelengths = np.linspace(_wavelengthRange[0],_wavelengthRange[1],nw)
    _qD = [AutoMieQ(_m,wavelength,diameter) for wavelength in wavelengths]
  elif type(_m) in [np.ndarray,list,tuple] and len(_wavelengthRange)==len(_m):
    wavelengths=_wavelengthRange
    _qD = [MieQ(emm,wavelength,diameter) for emm,wavelength in zip(_m,wavelengths)]
  else:
    warnings.warn("Error: the size of the input data is mismatched. Please examine your inputs and try again.")
    return None

  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return wavelengths, qext, qsca, qabs, g, qpr, qback, qratio

def MieQ_withSizeParameterRange(m, nMedium=1.0, xRange=(1,10), nx=1000, logX=False):
  """Calculate optical properties over a range of size parameters using Mie theory.

  Parameters
  ----------
  m : complex
      The complex refractive index of the particle material.
  nMedium : float, optional
      The real refractive index of the medium. Default is 1.0.
  xRange : tuple, optional
      A tuple specifying the range of size parameters (min, max). Default is (1, 10).
  nx : int, optional
      The number of size parameters to evaluate within the given range. Default is 1000.
  logX : bool, optional
      If `True`, the size parameters are spaced logarithmically between the minimum and maximum values. 
      If `False`, the spacing is linear. Default is `False`.

  Returns
  -------
  tuple
      Returns a tuple containing the following numpy arrays:
      - xValues: The size parameters evaluated.
      - qext: The extinction efficiency for each size parameter.
      - qsca: The scattering efficiency for each size parameter.
      - qabs: The absorption efficiency for each size parameter.
      - g: The asymmetry parameter for each size parameter.
      - qpr: The radiation pressure efficiency for each size parameter.
      - qback: The backscattering efficiency for each size parameter.
      - qratio: The ratio of backscattering to total scattering efficiency.

  Notes
  -----
  Given a range of size parameters (x) and the complex refractive index (m) of the scattering particle 
  relative to the medium, this function computes various optical properties using Mie theory.
  The size parameter "x" is defined as \\(x = 2 \\pi r / \\lambda\\), where "r" is the particle radius, 
  and "\\(\\lambda\\)" is the wavelength of the incident light in the medium.

  Examples
  --------
    ::

      import numpy as np
      m = 1.59 + 0.005j  # Complex refractive index of the particle
      nMedium = 1.33  # Refractive index of the medium
      xRange = (0.1, 100)  # Range of size parameters
      nx = 500  # Number of size parameters to evaluate
      logX = True  # Use logarithmic spacing

      # Compute the optical properties
      xValues, qext, qsca, qabs, g, qpr, qback, qratio = MieQ_withSizeParameterRange(m, nMedium, xRange, nx, logX)

      # Plot the extinction efficiency as a function of size parameter
      import matplotlib.pyplot as plt
      plt.loglog(xValues, qext)
      plt.xlabel('Size parameter (x)')
      plt.ylabel('Extinction efficiency (Qext)')
      plt.show()

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#MieQ_withSizeParameterRange
  nMedium = nMedium.real
  m /= nMedium
  _xRange = tuple([x*nMedium for x in xRange]) # I think
  if logX:
    xValues = list(np.logspace(np.log10(_xRange[0]),np.log10(_xRange[1]),nx))
  else:
    xValues = list(np.linspace(_xRange[0],_xRange[1], nx))
  dValues = [1000*x/np.pi for x in xValues]
  _qD = [AutoMieQ(m,1000,d) for d in dValues]
  qext = np.array([q[0] for q in _qD])
  qsca = np.array([q[1] for q in _qD])
  qabs = np.array([q[2] for q in _qD])
  g = np.array([q[3] for q in _qD])
  qpr = np.array([q[4] for q in _qD])
  qback = np.array([q[5] for q in _qD])
  qratio = np.array([q[6] for q in _qD])
  return xValues, qext, qsca, qabs, g, qpr, qback, qratio

def Mie_Lognormal(m,wavelength,geoStdDev,geoMean,numberOfParticles,nMedium=1.0, numberOfBins=10000,lower=1,upper=1000,gamma=[1],returnDistribution=False,decomposeMultimodal=False,asDict=False):
  """Perform Mie scattering calculations for lognormally distributed particles.

  Parameters
  ----------
  m : complex
      The refractive index of the particles.
  wavelength : float
      The wavelength of the incident light, in the same units as geoMean and lower/upper bounds.
  geoStdDev : float or sequence
      The geometric standard deviation(s) of the particles. Can be a single value or a list/tuple/numpy.ndarray for multimodal distributions.
  geoMean : float or sequence
      The geometric mean diameter(s) of the particles. Can be a single value or a list/tuple/numpy.ndarray for multimodal distributions.
  numberOfParticles : int
      The total number of particles.
  nMedium : float, optional
      The refractive index of the medium. Defaults to 1.0.
  numberOfBins : int, optional
      The number of bins for the discretization of the particle size distribution. Defaults to 10000.
  lower : int, optional
      The lower bound of particle diameters. Defaults to 1.
  upper : int, optional
      The upper bound of particle diameters. Defaults to 1000.
  gamma : list, optional
      The distribution of particles across different modes. Defaults to [1], implying a single mode.
  returnDistribution : bool, optional
      If True, the distribution data is also returned. Defaults to False.
  decomposeMultimodal : bool, optional
      If True and if the distribution is multimodal, return the contributions from each mode separately. Defaults to False.
  asDict : bool, optional
      If True, the optical properties are returned as a dictionary. Defaults to False.

  Returns
  -------
  tuple or dict
      Mie scattering properties (Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio) optionally along with the particle size distribution (dp, ndp) and, if multimodal and decomposeMultimodal is True, individual mode contributions (ndpi). The output is a tuple unless asDict is True, in which case it is a dictionary.

  Raises
  ------
  Warning
      If the distribution may not be compact on the specified interval or if not enough parameters are specified for each mode.

  Notes
  -----
  The function optionally warns the user if the specified particle size distribution may not be compact within the given interval or if there aren't enough parameters to fully specify each mode in the case of a multimodal distribution. This function is particularly useful in the analysis and simulation of light scattering by particles following a lognormal size distribution.

  """    
  #  http://pymiescatt.readthedocs.io/en/latest/forward.html#Mie_Lognormal
  nMedium = nMedium.real
  m /= nMedium
  wavelength /= nMedium
  ithPart = lambda gammai, dp, dpgi, sigmagi: (gammai/(np.sqrt(2*np.pi)*np.log(sigmagi)*dp))*np.exp(-(np.log(dp)-np.log(dpgi))**2/(2*np.log(sigmagi)**2))
  dp = np.logspace(np.log10(lower),np.log10(upper),numberOfBins)
  if all([type(x) in [list, tuple, np.ndarray] for x in [geoStdDev, geoMean]]):
    # multimodal
    if len(gamma)==1 and (len(geoStdDev)==len(geoMean)>1):
      # gamma is distributed equally among modes
      gamma = [1 for x in geoStdDev]
      gamma = [float(x/np.sum(gamma)) for x in gamma]
      ndpi = [numberOfParticles*ithPart(g,dp,dpg,sg) for g,dpg,sg in zip(gamma,geoMean,geoStdDev)]
      ndp = np.sum(ndpi,axis=0)
    elif len(gamma)==len(geoStdDev)==len(geoMean):
      # gamma is fully specified for each mode
      gamma = [float(x/np.sum(gamma)) for x in gamma]
      ndpi = [numberOfParticles*ithPart(g,dp,dpg,sg) for g,dpg,sg in zip(gamma,geoMean,geoStdDev)]
      ndp = np.sum(ndpi,axis=0)
    else:
      # user problem
      warnings.warn("Not enough parameters to fully specify each mode.")
      return None
  else:
    # unimodal
    decomposeMultimodal = False
    ndp = numberOfParticles*ithPart(1,dp,geoMean,geoStdDev)
  if ndp[-1]>np.max(ndp)/100 or ndp[0]>np.max(ndp)/100:
    warnings.warn("Warning: distribution may not be compact on the specified interval. Consider using a higher upper bound.")
  Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio = Mie_SD(m,wavelength,dp,ndp,SMPS=False)
  if returnDistribution:
    if decomposeMultimodal:
      if asDict==True:
        return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio), dp, ndp, ndpi
      else:
        return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio, dp, ndp, ndpi
    elif asDict==True:
      return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio), dp, ndp
    else:
      return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio, dp, ndp
  elif asDict==True:
    return dict(Bext=Bext, Bsca=Bsca, Babs=Babs, bigG=bigG, Bpr=Bpr, Bback=Bback, Bratio=Bratio)
  else:
    return Bext, Bsca, Babs, bigG, Bpr, Bback, Bratio
