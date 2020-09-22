import numpy as np
import ROOT as r
import math
import cmath
import os,sys
from scipy.integrate import quad, dblquad

from darkphoton import *

# proton mass
mProton = 0.938272081 # GeV/c - PDG2016
protonEnergy=400. # GeV/c

def protonMomentum(E=400.):
    return math.sqrt(E*E - mProton*mProton)

# useful functions
def energy(p,m):
    """ Compute energy from momentum and mass """
    return math.sqrt(p*p + m*m)

def rhoFormFactor(m):
    """ From https://arxiv.org/abs/0910.5589 """
    #f1prho = (3.2*0.77**2/5.03)/(0.77**2-m**2-0.77*0.15j)
    #f1pomega = (15.9*0.77**2/17.1)/(0.77**2-m**2-0.77*0.0085j)
    #return abs(f1prho+f1pomega)
    #simple Breit-Wigner
    #f1ra=math.sqrt(0.057412759)
    #f1rb=math.sqrt(0.30473311)
    #f1rc=math.sqrt(0.69790624)
    #f1wa=math.sqrt(0.00320847199)
    #f1wb=math.sqrt(0.30473311)
    #f1wc=math.sqrt(0.69790624)
    #code from Inar: https://github.com/pgdeniverville/BdNMC/blob/master/src/Proton_Brem_Distribution.cpp
    f1ra = 0.6165340033101271*0.77**2
    f1rb = 0.22320420111672623*1.25**2
    f1rc = -0.33973820442685326*1.45**2
    f1wa = 1.0117544786579074*0.77**2
    f1wb = -0.8816565944110686*1.25**2
    f1wc = 0.3699021157531611*1.45**2
    f1prho = f1ra/(0.77**2-m**2-0.77*0.15j)
    f1prhop = f1rb/(1.25**2-m**2-1.25*0.3j)
    f1prhopp = f1rc/(1.45**2-m**2-1.45*0.5j)
    f1pomega = f1wa/(0.77**2-m**2-0.77*0.0085j)
    f1pomegap = f1wb/(1.25**2-m**2-1.25*0.3j)
    f1pomegapp = f1wc/(1.45**2-m**2-1.45*0.5j)
    return abs(f1prho+f1prhop+f1prhopp+f1pomega+f1pomegap+f1pomegapp)

    #f1prho = (f1ra*0.77**2)**2/((m**2-0.77**2)**2+0.77**2*0.15**2)
    #f1prhop = (f1rb*1.25**2)**2/((m**2-1.25**2)**2+(1.25*0.3)**2)
    #f1prhopp = (f1rc*1.45**2)**2/((m**2-1.45**2)**2+(1.45*0.5)**2)
    #f1pomega = (f1wa*0.77**2)**2/((m**2-0.77**2)**2+(0.77*0.0085)**2)
    #f1pomegap = (f1wb*1.25**2)**2/((m**2-1.25**2)**2+(1.25*0.3)**2)
    #f1pomegapp = (f1wc*1.45**2)**2/((m**2-1.45**2)**2+(1.45*0.5)**2)
    #return math.sqrt(f1prho+f1prhop+f1prhopp+f1pomega+f1pomegap+f1pomegapp)
    #return math.sqrt(f1prho)+math.sqrt(f1pomega)
    #return f1prho+f1pomega

def penaltyFactor(m):
    """ Penalty factor for high masses - dipole form factor in the proton-A' vertex """
    """ m in GeV """
    if m*m>0.71:
        return math.pow(m*m/0.71,-4)
    else:
        return 1

def zeta(E, p, theta):
    """ Fraction of the proton momentum carried away by the paraphoton in the beam direction """
    return p / (protonMomentum(E) * math.sqrt(theta*theta + 1.))


def pTransverse(E, p, theta):
    """ Paraphoton transverse momentum in the lab frame """
    return protonMomentum(E)*theta*zeta(E,p,theta)


def ptSquare(E, p, theta):
    """ Square paraphoton transverse momentum in the lab frame """
    return pow(pTransverse(E,p,theta), 2.)


def H(E, p, theta, mDarkPhoton):
    """ A kinematic term """
    return ptSquare(E,p,theta) + (1.-zeta(E,p,theta))*mDarkPhoton*mDarkPhoton + pow(zeta(E,p,theta),2.)*mProton*mProton


def wba(E, p, theta, mDarkPhoton, epsilon):
    """ Cross section weighting function in the Fermi-Weizsaeker-Williams approximation """
    const = epsilon*epsilon*alphaQED / (2.*math.pi*H(E,p,theta,mDarkPhoton))

    h2 = pow(H(E,p,theta,mDarkPhoton),2.)
    oneMinusZSquare = pow(1.-zeta(E,p,theta),2.)
    mp2 = mProton*mProton
    mA2 = mDarkPhoton*mDarkPhoton

    p1 = (1. + oneMinusZSquare) / zeta(E,p,theta)
    p2 = ( 2. * zeta(E,p,theta) * (1.-zeta(E,p,theta)) * ( (2.*mp2 + mA2)/ H(E,p,theta,mDarkPhoton) 
            - pow(zeta(E,p,theta),2.)*2.*mp2*mp2/h2 ) )
    #p3 = 2.*zeta(p,theta)*(1.-zeta(p,theta))*(zeta(p,theta)+oneMinusZSquare)*mp2*mA2/h2
    p3 = 2.*zeta(E,p,theta)*(1.-zeta(E,p,theta))*(1+oneMinusZSquare)*mp2*mA2/h2
    p4 = 2.*zeta(E,p,theta)*oneMinusZSquare*mA2*mA2/h2
    return const*(p1-p2+p3+p4)


def sigma(s): # s in GeV^2 ---> sigma in mb
    """ Parametrisation of sigma(s) """
    a1 = 35.45
    a2 = 0.308
    a3 = 28.94
    a4 = 33.34
    a5 = 0.545
    a6 = 0.458
    a7 = 42.53
    p1 = a2*pow(math.log(s/a3),2.) 
    p2 = a4*pow((1./s),a5)
    p3 = a7*pow((1./s),a6)
    return a1 + p1 - p2 + p3


def es(E, p, mDarkPhoton):
    """ s(p,mA) """
    return 2.*mProton*(energy(protonMomentum(E),mProton)-energy(p,mDarkPhoton))


def sigmaRatio(E, p, mDarkPhoton):
    """ sigma(s') / sigma(s) """
    return sigma(es(E,p,mDarkPhoton)) / sigma(2.*mProton*energy(protonMomentum(E),mProton))


def dNdZdPtSquare(E, p, mDarkPhoton, theta, epsilon):
    """ Differential A' rate per p.o.t. as a function of Z and Pt^2 """
    return sigmaRatio(E,p,mDarkPhoton)*wba(E,p,theta,mDarkPhoton,epsilon)


def dPt2dTheta(E,p, theta):
    """ Jacobian Pt^2->theta """
    z2 = pow(zeta(E,p,theta),2.)
    return 2.*theta*z2*protonMomentum(E)*protonMomentum(E)


def dZdP(E, p, theta):
    """ Jacobian z->p """
    return 1./( protonMomentum(E)* math.sqrt(theta*theta+1.) )


def dNdPdTheta(p, theta, E, mDarkPhoton, epsilon):
    """ Differential A' rate per p.o.t. as a function of P and theta """
    diffRate = dNdZdPtSquare(E,p,mDarkPhoton,theta,epsilon) * dPt2dTheta(E,p,theta) * dZdP(E,p,theta)
    return math.fabs(diffRate) # integrating in (-pi, pi)...


def pMin(E,mDarkPhoton):
    return max(0.14*protonMomentum(E), mDarkPhoton)

def pMinFix(E,mDarkPhoton):
    return max(0.1*protonMomentum(E), mDarkPhoton)


def pMax(E,mDarkPhoton):
    #return min(0.86*protonMomentum, math.sqrt( (energy(protonMomentum,mProton)**2. - mDarkPhoton**2.) - mDarkPhoton**2.))
    return math.sqrt( (energy(protonMomentum(E),mProton)**2. - mDarkPhoton**2.) - mDarkPhoton**2.)

def pMaxFix(E,mDarkPhoton):
    return min(0.9*protonMomentum(E), math.sqrt( (energy(protonMomentum(E),mProton)**2. - mDarkPhoton**2.) - mDarkPhoton**2.))


def prodRate(E, mDarkPhoton, epsilon, tmin = -0.5 * math.pi, tmax = 0.5 * math.pi, pmin = -1, pmax = -1):
    if pmin==-1: pmin = pMin(E,mDarkPhoton)
    if pmax==-1: pmax = pMax(E,mDarkPhoton)
    """ dNdPdTheta integrated over p and theta """
    integral = dblquad( dNdPdTheta, # integrand
                        tmin, tmax, # theta boundaries (2nd argument of integrand)
#                        lambda x: pMin(E,mDarkPhoton), lambda x: pMax(E,mDarkPhoton), # p boundaries (1st argument of integrand)
                        lambda x: pmin, lambda x: pmax, # p boundaries (1st argument of integrand)
                        args=(E,mDarkPhoton, epsilon) ) # extra parameters to pass to integrand
    return integral[0]

# total production rate of A'
#norm = prodRate(1.1,3.e-7) #mDarkPhoton,epsilon)
# number of A' produced
# numDarkPhotons = int(math.floor(norm*protonFlux))
#
# print
# print "Epsilon \t %s"%epsilon
# print "mDarkPhoton \t %s"%mDarkPhoton
#print "A' production rate per p.o.t: \t %.8g"%norm
# print "Number of A' produced in SHiP: \t %.8g"%numDarkPhotons


def normalisedProductionPDF(E, p, theta, mDarkPhoton, epsilon, norm):
    """ Probability density function for A' production in SHIP """
    return (1. / norm) * dNdPdTheta(p, theta, E, mDarkPhoton, epsilon)


def hProdPDF(E,mDarkPhoton, epsilon, norm, binsp, binstheta, tmin = -0.5 * math.pi, tmax = 0.5 * math.pi, pmin = -1, pmax = -1, suffix=""):
    if pmin==-1: pmin = pMin(E,mDarkPhoton)
    if pmax==-1: pmax = pMax(E,mDarkPhoton)
    """ Histogram of the PDF for A' production in SHIP """
    angles = np.linspace(tmin,tmax,binstheta).tolist()
    anglestep = 2.*(tmax - tmin)/binstheta
    momentumStep = (pmax-pmin)/(binsp-1)
    momenta = np.linspace(pmin,pmax,binsp,endpoint=False).tolist()
    hPDF = r.TH2F("hPDF_p%s_eps%s_m%s"%(E,epsilon,mDarkPhoton) ,"hPDF_p%s_eps%s_m%s"%(E,epsilon,mDarkPhoton),
        binsp,pmin-0.5*momentumStep,pmax-0.5*momentumStep,
        binstheta,tmin-0.5*anglestep,tmax-0.5*anglestep)
    hPDF.SetTitle("PDF for A' production (Ep=%s GeV, m_{A'}=%s GeV, #epsilon =%s)"%(E,mDarkPhoton,epsilon))
    hPDF.GetXaxis().SetTitle("P_{A'} [GeV]")
    hPDF.GetYaxis().SetTitle("#theta_{A'} [rad]")
    hPDFtheta = r.TH1F("hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
        "hPDFtheta_eps%s_m%s"%(epsilon,mDarkPhoton),
        binstheta,tmin-0.5*anglestep,tmax-0.5*anglestep)
    hPDFp = r.TH1F("hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
        "hPDFp_eps%s_m%s"%(epsilon,mDarkPhoton),
        binsp,pmin-0.5*momentumStep,pmax-0.5*momentumStep)
    hPDFp.GetXaxis().SetTitle("P_{A'} [GeV]")
    hPDFtheta.GetXaxis().SetTitle("#theta_{A'} [rad]")
    for theta in angles:
        for p in momenta:
            w = normalisedProductionPDF(E,p,theta,mDarkPhoton,epsilon,norm)
            hPDF.Fill(p,theta,w)
            hPDFtheta.Fill(theta,w)
            hPDFp.Fill(p,w)
    hPdfFilename = sys.modules['__main__'].outputDir+"/ParaPhoton_eps%s_m%s%s.root"%(epsilon,mDarkPhoton,suffix)
    outfile = r.TFile(hPdfFilename,"recreate")
    #weight = hPDF.Integral("width")
    #print "Weight = %3.3f"%weight
    #hPDF.Scale(1./weight)
    hPDF.Write()
    hPDFp.Write()
    hPDFtheta.Write()
    outfile.Close()
    del angles
    del momenta
    return hPDF


