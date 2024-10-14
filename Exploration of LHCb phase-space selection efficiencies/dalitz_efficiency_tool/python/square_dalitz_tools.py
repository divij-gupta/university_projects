import numpy as np
from ROOT import TLorentzVector

def get_q(e1, m1):
    """
    calculate q (momentum) in particle rest frame.

    :param e1: energy of particle  m1
    :param m1: invariant mass of particle m1
    :returns: q
    """

    return np.sqrt(np.power(e1, 2) - np.power(m1, 2))

def get_theta(m12, m13, m1=0.497614, m2=0.13957061, m3=0.1395706, m123=1.86483):
    """
    calculate theta (angle between d1 d3 in the d12 rest frame)

    :param m12: mass of m12
    :param m13: mass of m13
    :param m1: mass of m1. defaults to Ks
    :param m2: mass of m2. defaults to pion
    :param m3: mass of m3. defaults to pion
    :param m123: mass of m123 defaults to D0
    :returns: h_12
    """

    #calculate invariant masses from masses
    s123 = np.power(m123, 2)
    s12 = np.power(m12, 2)
    s13 = np.power(m13, 2)
    s1 = np.power(m1, 2)
    s2 = np.power(m2, 2)
    s3 = np.power(m3, 2)

    #centre of mass energies for d1 and d3 in m12 frame
    e_1 = (s12 - s2 + s1)/(2*m12)
    e_3 = (s123 - s12 - s3)/(2*m12)

    #check energy is physical
    zeros_1 = np.where(e_1 < m1)[0]
    zeros_3 = np.where(e_3 < m3)[0]
    zeros = np.append(zeros_1, zeros_3)

    #calculate momentum of d1 and d3
    q_1 = get_q(e_1, m1)
    q_3 = get_q(e_3, m3)

    #cosine of helicity angle
    cos_h_12 = -(s13 - s1 - s3 - 2.0*e_1*e_3)/(2.0*q_1*q_3)

    #check physical boundaries
    cos_h_12 = np.where(cos_h_12 > 1, 1, cos_h_12)
    cos_h_12 = np.where(cos_h_12 < -1, -1, cos_h_12)

    h_12 = np.arccos(cos_h_12)
    h_12[zeros] = 0

    return h_12

def get_thetaprime(theta):
    """
    calculate theta' for square dalitz formalism.

    :param theta: theta
    :returns: theta'
    """

    return (1/np.pi)*theta

def get_kinlimits(m1=0.497614, m2=0.13957061, m3=0.1395706, m123=1.86483):
    """
    return kinematic limits for the dalitz space. defaults to the s13 s23 space. Specific to KSpipi

    :param m1: mass of m1. defaults to Ks
    :param m2: mass of m2. defaults to pion
    :param m3: mass of m3. defaults to pion
    :param m123: mass of m123 defaults to D0
    :returns: kinimatic limits Mmin mMax
    """

    m12_min = np.abs(m1 + m2)
    m12_max = np.abs(m123 - m3)

    return m12_min, m12_max

def get_mprime(m12, m1=0.497614, m2=0.13957061, m3=0.1395706, m123=1.86483):
    """
    calculate m' for square dalitz formalism.

    :param m12: mass m12
    :param m1: mass of m1. defaults to Ks
    :param m2: mass of m2. defaults to pion
    :param m3: mass of m3. defaults to pion
    :param m123: mass of m123 defaults to D0
    :returns: m'
    """

    m12_min, m12_max = get_kinlimits(m1, m2, m3, m123)
    a = (m12 - m12_min)/(m12_max - m12_min)                               #reduce clutter on the return

    return (1/np.pi)*np.arccos((2*a) - 1)

def get_jacobian(m12, mprime, thetaprime, m1=0.497614, m2=0.13957061, m3=0.1395706, m123=1.86483):

    m12_min, m12_max = get_kinlimits(m1, m2, m3, m123)

    a = -(np.pi/2) * np.sin(np.pi * mprime) * (m12_max - m12_min)
    b = -np.pi * np.sin(np.pi * thetaprime)

    # calculate invariant masses from masses
    s12 = np.power(m12, 2)
    s123 = np.power(m123, 2)
    s1 = np.power(m1, 2)
    s2 = np.power(m2, 2)
    s3 = np.power(m3, 2)

    # centre of mass energies for d1 and d3 in m12 frame
    e_1 = (s12 - s2 + s1) / (2 * m12)
    e_3 = (s123 - s12 - s3) / (2 * m12)

    p = get_q(e_3, m3)                                                   # mom m3
    q = get_q(e_1, m1)                                                   # mom m1

    e_123 = np.sqrt(s123 + p**2)
    gamma = e_123 / m123
    beta = p / e_123

    q_boost = np.sqrt((beta * gamma * e_1)**2 + q**2)
    p_boost = gamma * p - gamma * beta * e_3

    # return 4 * p_boost * q_boost * m12 * a * b
    return 4 * p * q * m12 * a * b

def make_vector(df, prefix, particle):
    """
    make TLorentzVectors for a given particle

    :param df: input dataframe containing PX PY PZ and E for a given particle
    :prefix: decaytree fitter prefix
    :particle: particle name (h1, KS etc)
    :returns: dataframe with TLorentzVector column for requested particle
    """

    name = f"TLorentzVector_{particle}"
    if name in df.columns:
        print(f"...{name} already exists!")
    else:
        print(f"...generating {name}")
        df[name] = df.apply(lambda row : TLorentzVector(row[f"{prefix}_{particle}_PX"],
                                                        row[f"{prefix}_{particle}_PY"],
                                                        row[f"{prefix}_{particle}_PZ"],
                                                        row[f"{prefix}_{particle}_E"]), axis = 1)
        print("...done!")

def make_dalitz(df, h1, h2, h3):
    """
    make square dalitz variables with respect to h1 h2. Converting to MeV! 

    :param df: input dataframe containing TLorentzVectors for h1-3
    :h1: name of h1 particle
    :h2: name of h2 particle
    :h3: name of h3 particle
    :returns: dataframe with square dalitz variables
    """

    df["mpipi"] = df.apply(lambda row : get_mprime((row[f"TLorentzVector_{h1}"]
                                                   +row[f"TLorentzVector_{h2}"]).M()/1000),
                            axis = 1)
    df["costh"] = df.apply(lambda row : get_thetaprime(get_theta((row[f"TLorentzVector_{h1}"]
                                                                  +row[f"TLorentzVector_{h2}"]).M()/1000,
                                                                (row[f"TLorentzVector_{h1}"]
                                                                  +row[f"TLorentzVector_{h3}"]).M()/1000)),
                           axis = 1)
