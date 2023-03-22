#!/usr/bin/env python

'''

GOAL:
- Measure completeness corrections of LDP spectroscopic catalogs.


USAGE:
- in ipython


REQUIRED MODULES:


**************************
written by Gregory Rudnick
Deceber 2022
**************************

'''

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import constants as c 
from astropy import units as u
from astropy.table import Table, Column
import os
import sys
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
from matplotlib.image import NonUniformImage
#import warnings
#warnings.filterwarnings("ignore")

homedir = os.getenv("HOME")
ldppath = homedir + '/EDisCS/Catalogs/LDP/'

'''
methods to include

TO DO

check how CMD histograms compare from field to field

write out global completeness and cluster-by-cluster completeness.  Do
this only for the CMD version


------------------------
**radial completeness

input: total cat, sample flag
cluster: default all, can input name

Make radial histogram for all photometric sources

make radial histogram for all targets with a Q=3 or 4 spectrum (or real spectrum)

take ratio of these to compute completeness

returns
completeness radial 1D array with bin edges

---------------------------
**assign completeness weights

input: total cat, sample flag, CMD completeness, radial completeness

for every galaxy with a redshift find corresponding CMD bin and assign completeness 

for every galaxy with a redshift find corresponding radial bin and assign completeness 

returns
radial and CMD completeness for every source

'''
class complete:

    
    def __init__(self, ldppath):
        catname = ldppath + 'v7.2/megacat_v7.2a.fits'
        self.catdat, self.cathdr = fits.getdata(catname,header = True)
        self.cat = Table(self.catdat)
        #hdul = fits.open(catname)
        #cat = hdul[1].data
        #cathdr = hdul[0].header
        

    def mwextcor(self):

        '''

        PURPOSE: 

        correct the full flux catalog for galactic extinction
        #use the following formula F_obs(Lambda) = F_0(Lambda) * 10^(-0.4 * EBV * K(Lambda))
        #Calzetti+97 https://arxiv.org/pdf/astro-ph/9706121.pdf

        ebv values have already been tabulated for each source

        INPUT:
        The full flux catalog

        OUTPUT:
        A new astropy table that contains corrected fluxs with the '_c' extension for every aperture.


        '''
        self.klam = {'B' : 4.125,
                    'V' : 3.060,
                    'R' : 2.548,
                    'I' : 1.993,
                    'z' : 1.384,
                    'K' : 0.367}

        #the different apertures
        apertype = ['1', '2','3','iso','auto']
        apertypek = ['1', '2', '3']

        #the different variants of the K-band
        ktype = ['','_UKIRT']

        #initialize new astropy Table with new columns 
        self.nc = Table()
        
        #loop through every band
        for band in self.klam:

            #loop through both variants of the K-band
            if band == 'K':
                #loop through apertypek as there are different
                #aperture types for the k-band (not sure why but there
                #are)
                for aper in apertypek:
                    #loop through NEWFIRM ('') and UKIRT ('_UKIRT') extensions
                    for ext in ktype:
                        self.fluxcor(band,aper,ext)
            else:
                #loop through apertype
                ext = ''
                for aper in apertype:
                    self.fluxcor(band,aper,ext)

                
    def fluxcor(self,band,aper,ext):
        '''
        correct the fluxes and errors for galactic extinction for a
        single band and aperture combination.

        The 'ext' parameter is only important for the UKIRT data

        '''

        #make a new variable that contains the corrected flux
        fluxstr = 'f' + band + aper + ext
        fluxstrebvcor = fluxstr + '_c'
        self.nc[fluxstrebvcor] = self.cat[fluxstr] * 10**(0.4 * self.cat['ebv'] * self.klam[band])

        #make a new variable that contains the corrected flux errors
        errfluxstr = 'f' + band + aper + 'err' + ext
        errfluxstrebvcor = 'f' + band + aper + 'err' + ext + '_c'
        self.nc[errfluxstrebvcor] = self.cat[errfluxstr] * 10**(0.4 * self.cat['ebv'] * self.klam[band])

        #replace values with no data with their original values
        noflux = np.where(self.cat[errfluxstr] < 0.)
        self.nc[fluxstrebvcor][noflux] = self.cat[fluxstr][noflux]
        self.nc[errfluxstrebvcor][noflux] = self.cat[errfluxstr][noflux]

    def galsel(self):

        '''

        select the subset of galaxies with good photometry and
        good redshifts by creating flags apply quality cuts to make
        sampleflag

        RETURNS
        sampleflag
        '''
            
        #self.mRAUTO = -2.5 * np.log10(self.cat.fRauto) + 23.9
        self.mRAUTO = -2.5 * np.log10(self.nc['fRauto_c']) + 23.9
        self.starflag = (self.cat['class_StarR']>0.98)
        #select galaxies with high quality photometry 
        self.photflag = (self.cat['wK_both']>0.3) & (self.cat['sexflagB']<=3) & (self.cat['sexflagV']<=3) &  (self.cat['sexflagR']<=3)&  (self.cat['sexflagI']<=3) & ((self.cat['sexflagK'] <=3) | (self.cat['sexflagK_UKIRT'] <=3)) & (~self.starflag)
        print('nphot_all = ', np.count_nonzero(self.photflag))

        #select galaxies with high quality spectroscopic redshifts 
        self.zflag = ((self.cat['zGAL_new_ORIGIN']=='LDP') & (self.cat['Q'] >=3)) | (self.cat['zGAL_new_ORIGIN']=='FORS') | (self.cat['zGAL_new_ORIGIN']=='HS')
        print('nspec_all = ', np.count_nonzero(self.zflag))
        
        #add spectroscopic completness limit
        self.goodphotflag = self.photflag & (self.mRAUTO < 22.9)
        self.goodspecflag = self.goodphotflag & self.zflag
        #self.goodspecflag = self.zflag  & (self.mRAUTO < 22.9)
        print('nphot = ',np.count_nonzero(self.goodphotflag))
        print('nspec =',np.count_nonzero(self.goodspecflag))

    def redcheck(self):
        '''
        check redenning correction


        '''

        V_R_c = -2.5 * np.log10(self.nc['fV1_c'] / self.nc['fR1_c'])
        V_R = -2.5 * np.log10(self.cat['fV1'] / self.cat['fR1'])

        mRAUTO = -2.5 * np.log10(self.cat['fRauto']) + 23.9
        mRAUTOc = -2.5 * np.log10(self.nc['fRauto_c']) + 23.9

        fig,(ax1,ax2) = plt.subplots(1,2,figsize = (8,4))  

        ax1.scatter(V_R_c[self.goodphotflag], V_R_c[self.goodphotflag] - V_R[self.goodphotflag], alpha=0.1)
        ax1.set_xlim(-1,2)
        ax1.set_ylabel(r'(V-R)$_c$ - (V-R)')
        ax1.set_xlabel('(V-R)_c')

        ax2.scatter(mRAUTO[self.goodphotflag], mRAUTOc[self.goodphotflag] - mRAUTO[self.goodphotflag],alpha=0.1)
        ax2.set_ylabel(r'RAUTO$_c$ - RAUTO')
        ax2.set_xlabel('RAUTO')

        fig,ax = plt.subplots(1,1,figsize = (4,4))
        ax.scatter(mRAUTOc[self.goodphotflag],V_R_c[self.goodphotflag] - V_R[self.goodphotflag], alpha=0.1)
        ax.set_ylabel(r'(V-R)$_c$ - (V-R)')
        ax.set_xlabel('RAUTO')

    def CMD_compl(self, cluster='all', debug = False):
        '''CMD completeness

        INPUT: 
        total cat
        sample flag
        cluster: default all, can input name

        WHAT IT DOES:
        make 2D binned histogram of V-R vs. R for all photometric sources

        make 2D binned histogram of V-R vs. R for all targets with a Q=3 or 4 spectrum (or real spectrum)

        take ratio of these to compute completeness

        RETURNS
        completeness CMD 2D array with bin edges
        '''

        #select the right cluster
        if(cluster=='all'):
            self.clustflag = (self.cat['FIELD']!='None')
        else:
            self.clustflag = (self.cat['FIELD']==cluster)

        print('numclust = ', np.count_nonzero(self.clustflag))
        #print('numphot = ', np.count_nonzero(self.goodphotflag))
        #print('numspec = ', np.count_nonzero(self.goodspecflag))

        self.goodphot_plotflag = self.goodphotflag & self.clustflag
        self.goodspec_plotflag = self.goodspecflag & self.clustflag
        print('numphotclust = ', np.count_nonzero(self.goodphot_plotflag))
        print('numspecclust = ', np.count_nonzero(self.goodspec_plotflag))

        #make plot of V-R vs. R
        #self.V_R = -2.5 * np.log10(self.cat.fV1 / self.cat.fR1)
        self.V_R = -2.5 * np.log10(self.nc['fV1_c'] / self.nc['fR1_c'])

        fig,(ax1,ax2) = plt.subplots(1,2,figsize = (8,4))  

        #plt.figure(figsize=(8,6))
        #ax=plt.gca()
        #fig, ax = plt.subplots(nrows=1, ncols=2)
        ax1.scatter(self.mRAUTO[self.goodphot_plotflag], self.V_R[self.goodphot_plotflag],alpha=0.2)
        #plt.scatter(self.mRAUTO[self.goodphotflag], self.V_R[self.goodphotflag])

        ax1.set_xlim([23.3,15])
        ax1.set_ylim([-1,2])
    
        ax1.set_ylabel(r'V-R',fontsize=20)
        ax1.set_xlabel(r'$R_{AUTO}$',fontsize=20)
        ax1.set_title('Phot',fontsize=20)
        ax1.text(18,1.6,s=cluster,fontsize=20)
    
        #ax1.legend(loc=2,fontsize=18)
        ax1.tick_params(axis='both', which='major', labelsize=15)

        ax2.scatter(self.mRAUTO[self.goodspec_plotflag], self.V_R[self.goodspec_plotflag],alpha=0.5)
        #plt.scatter(self.mRAUTO[self.goodphotflag], self.V_R[self.goodphotflag])

        ax2.set_xlim([23.3,15])
        ax2.set_ylim([-1,2])
    
        #ax2.set_ylabel(r'V-R',fontsize=20)
        ax2.set_xlabel(r'$R_{AUTO}$',fontsize=20)
        ax2.set_title('Spec-z',fontsize=20)

        #ax2.legend(loc=2,fontsize=18)

        fig.subplots_adjust(wspace = 0.15)

        ax2.tick_params(axis='both', which='major', labelsize=15)

        figname = '../Plots/cmd_' + cluster + '.png'
        fig.savefig(figname)
        

        ##############
        #compute 2D histogram 
        #set up panel for histograms
        fig,(ax1,ax2) = plt.subplots(1,2,figsize = (8,4))

        #set boundaries of histogram for completeness.  mRAUTO=22.9 is
        #the spectroscopic completeness as estimated in Just+19
        self.magedges = np.arange(15.9,23.9,0.5)
        self.coledges = np.arange(-0.5,2.0,0.25)

        #do it for photometric sources
        histphot, xedges,yedges = np.histogram2d(self.mRAUTO[self.goodphot_plotflag], self.V_R[self.goodphot_plotflag], bins=[self.magedges,self.coledges])
        histphot = histphot.T

        xphot, yphot = np.meshgrid(self.magedges,self.coledges)
        im1 = ax1.pcolormesh(xphot, yphot, np.log10(histphot))
        cb1 = fig.colorbar(im1,ax=ax1)
        cb1.set_label(label=r'log$_{10}$ N',fontsize=15)
        cb1.ax.tick_params(labelsize=15)
        ax1.set_xlim([23.3,15])
        ax1.set_ylim([-1,2])
    
        ax1.set_ylabel(r'V-R',fontsize=20)
        ax1.set_xlabel(r'$R_{AUTO}$',fontsize=20)
        ax1.set_title('Phot',fontsize=20)
        ax1.text(18,1.6,s=cluster,fontsize=20)

        ax1.tick_params(axis='both', which='major', labelsize=15)

        #now for spectroscopic sources
        histspec, xedges,yedges = np.histogram2d(self.mRAUTO[self.goodspec_plotflag], self.V_R[self.goodspec_plotflag], bins=[self.magedges,self.coledges])
        histspec = histspec.T

        xspec, yspec = np.meshgrid(self.magedges,self.coledges)
        im2 = ax2.pcolormesh(xspec, yspec, np.log10(histspec))
        cb2 = fig.colorbar(im2,ax=ax2)
        cb2.set_label(label=r'log$_{10}$ N',fontsize=15)
        cb2.ax.tick_params(labelsize=15)
        ax2.set_xlim([23.3,15])
        ax2.set_ylim([-1,2])
    
        #ax2.set_ylabel(r'V-R',fontsize=20)
        ax2.set_xlabel(r'$R_{AUTO}$',fontsize=20)
        ax2.set_title('Spec-z',fontsize=20)

        ax2.tick_params(axis='both', which='major', labelsize=15)

        figname = '../Plots/cmdhist_' + cluster + '.png'
        fig.savefig(figname)
        #plt.show()

        ###############################3
        #now create the completeness array
        self.completehist = histspec / histphot
        #find cells with few photometric sources
        zerophot = np.where(histphot<10)
        self.completehist[zerophot] = 0

        if debug:
            self.completehist[3,4] = 0.7
        
        fig,ax = plt.subplots(1,1,figsize = (4,4))
        im = ax.pcolormesh(xspec, yspec, self.completehist,vmin=0,vmax = 0.6)
        cb = fig.colorbar(im,ax=ax)
        #cb.set_label(label=r'log$_{10}$ N',fontsize=15)
        cb.ax.tick_params(labelsize=15)
        cb.set_label(label=r'$N_{spec}/N_{phot}$',fontsize=15)

        ax.set_xlim([23.3,15.9])
        ax.set_ylim([-0.5,1.75])
    
        ax.set_ylabel(r'V-R',fontsize=20)
        ax.set_xlabel(r'$R_{AUTO}$',fontsize=20)
        ax.set_title('Completeness',fontsize=20)
        ax.text(18,1.6,s=cluster,fontsize=20,color='white')

        ax.tick_params(axis='both', which='major', labelsize=15)

        figname = '../Plots/compl_' + cluster + '.png'
        fig.savefig(figname)
        plt.show()


    def CMD_compl_assign(self):
        '''CMD_compl_assign

        for every galaxy find its closes cell in the global and
        field-by-field completeness histogram and add the completeness
        value for that galaxy.

        OUTPUT:

        The original catalog with a new completeness column

        PROCESS:

        Loop through every CMD bin.  Find all galaxies in that bin and
        assign them the appropriate completeness
        '''

        #initialize weight column
        self.w_cmd = np.zeros(len(self.cat))
        #print(self.w_cmd)

        #loop over magnitude bins
        for imag,magbright in enumerate(self.magedges):
            #only go to second to last edge so that we can access the final bin and not beyond
            if imag < len(self.magedges) - 1:
                magfaint = self.magedges[imag+1]

                #loop over color bins
                for icol, colblue in enumerate(self.coledges):
                    if icol < len(self.coledges) - 1:
                        colred = self.coledges[icol + 1]

                        #find all galaxies in this 
                        #galbinflag = (self.mRAUTO <= magbright) & (self.mRAUTO > magfaint) & (self.V_R  >= colblue) & (self.V_R < colred)
                        galbinflag = (self.mRAUTO >= magbright) & (self.mRAUTO < magfaint) & (self.V_R  >= colblue) & (self.V_R < colred)
                        
                        #print('Ngal selected = ',np.count_nonzero(galbinflag))
                        #now assign the appropriate weight to every galaxy in that bin 
                        self.w_cmd[galbinflag] = self.completehist[icol,imag]
                        #print('mag = ',self.magedges[imag], '; magfaint = ', magfaint, 'magbright = ', magbright, '; col = ',self.coledges[icol], self.completehist[icol,imag])
                        #print('w_cmd[galbinflag] = ', self.w_cmd[galbinflag])

        self.cat['w_cmd'] = self.w_cmd

        plt.figure(figsize = (6,6))
        plt.scatter(self.mRAUTO[self.goodphot_plotflag], self.V_R[self.goodphot_plotflag],alpha=1, c=self.cat['w_cmd'][self.goodphot_plotflag])
        cb = plt.colorbar()
        cb.set_label(label=r'$w_{cmd}$',fontsize=15)

        #plt.scatter(self.mRAUTO[self.goodphotflag], self.V_R[self.goodphotflag])

        plt.xlim(23.3,15)
        plt.ylim(-1,2)
    
        plt.ylabel(r'V-R',fontsize=20)
        plt.xlabel(r'$R_{AUTO}$',fontsize=20)
        plt.title('Phot',fontsize=20)
        #plt.text(18,1.6,s=cluster,fontsize=20)
    
        #ax1.legend(loc=2,fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=15)
        figname = '../Plots/cmd_weights.png'
        plt.savefig(figname)
        

    def compl_write(self):
        '''compl_write

        write an updated catalog with completeness weights

        OUTPUT:

        A fits file with the added completeness column

        '''
        catcomplname = ldppath + 'v7.2/megacat_v7.2a_compl.fits'
        self.cat.write(catcomplname, overwrite=True)
        
        
if __name__ == '__main__':
    c = complete(ldppath)
    c.mwextcor()
    c.galsel()
