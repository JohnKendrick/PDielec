#!/usr/bin/python
from __future__ import print_function, division
import string
import re
import numpy as np
import math 
import os, sys

class Plotter:
    """Do some plotting for the dielectric functions"""
    def __init__(self):
        self.methods = []
        self.volume_fractions = []
        self.shapes = []
        self.shape_data = []
        self.plot_numbers = []
        self.traces = []
        self.absorption_coefficients = []
        self.molar_absorption_coefficients = []
        self.frequencies = []
        return

    def addDielectric(self, nplot, method, vf_type, shape, data, v, trace, absorption, molar_absorption) :
        """Add the trace of a dielectric to be plotted, add the absorption and the molar absorption"""
        if not nplot in self.plot_numbers:
            self.methods.append(method)
            self.volume_fractions.append(vf_type)
            self.shapes.append(shape)
            self.shape_data.append(data)
            self.plot_numbers.append(nplot)
            self.traces.append([])
            self.absorption_coefficients.append([])
            self.molar_absorption_coefficients.append([])
            self.frequencies.append([])
        self.frequencies[nplot].append(v)
        self.traces[nplot].append(trace)
        self.absorption_coefficients[nplot].append(absorption)
        self.molar_absorption_coefficients[nplot].append(molar_absorption)
        return 

    def printout(self, fd_csvfile) :
        """Print data to a csv file"""
        if fd_csvfile == 0:
            return
        # Print the header
        output = ",,,,,cm-1,thz"
        for method,shape,data,vf_type in zip(self.methods,self.shapes,self.shape_data,self.volume_fractions):
             strdata = str(data)
             strdata=strdata.replace(" ","")
             strdata=strdata.replace(","," ")
             strdata=strdata.replace("'","")
             strdata=strdata.replace('"',"")
             output += ","+method+"_"+shape+strdata+"_real_dielectric"+"("+vf_type+")"
             output += ","+method+"_"+shape+strdata+"_imag_dielectric"+"("+vf_type+")"
             output += ","+method+"_"+shape+strdata+"_absorption_coeficient_cm-1"+"("+vf_type+")"
             output += ","+method+"_"+shape+strdata+"_molar_absorption_coeficient_cm-1.L.mole-1"+"("+vf_type+")"
            # end loop over shape
        #end loop over method
        print(output, file=fd_csvfile)
        # Use the first frequency list
        for iv,v in enumerate(self.frequencies[0]):
            thz = v * 0.0299792458
            output = ",,,,,{:12f},{:12f}".format(v,thz)
            for dri,absorption,molar_absorption in zip(self.traces,self.absorption_coefficients,self.molar_absorption_coefficients):
                    dr = float(np.real(dri[iv]))
                    di = float(np.imag(dri[iv]))
                    K1 = absorption[iv]
                    K2 = molar_absorption[iv]
                    output += "," + ",".join("{:20.8f}".format(p) for p in [dr, di, K1, K2] )
                # end loop over shape
            # end loop over method
            print(output, file=fd_csvfile)
        # end loop over frequency
        return

    def plot(self,types) :
        for plot in types:
            if plot == 'imaginary':
                self.plotImaginary()
            elif plot == 'real':
                self.plotReal()
            elif plot == 'absorption':
                self.plotAbsorption()
            elif plot == 'molar_absorption':
                self.plotMolarAbsorption()
            elif plot == 'molarAbsorption':
                self.plotMolarAbsorption()
            elif plot == 'molarabsorption':
                self.plotMolarAbsorption()
            else:
                print('Unkown plot type {}'.format(plot))
                exit(1)
            # endif
        #endfor
        return

    def plotMolarAbsorption(self) :
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(self.molar_absorption_coefficients[-1]) 
        # loop over all methods except the last, TO
        for molar_absorption in self.molar_absorption_coefficients[:-1]:
            maxd = max(maxd, np.amax(molar_absorption) )
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        molar_absorption = np.array(self.molar_absorption_coefficients[-1])
        self.molar_absorption_coefficients[-1] = to_scale * molar_absorption
        plot_list = []
        plot_labels = []
        for method,shape,data,frequencies,molar_absorption,vf_type in zip(self.methods,self.shapes,self.shape_data,self.frequencies,self.molar_absorption_coefficients,self.volume_fractions):
           label = method+"_"+shape+str(data)+"("+vf_type+")"
           label=label.replace(" ","")
           label=label.replace(","," ")
           label=label.replace("'","")
           label=label.replace('"','')
           plot, = pl.plot(frequencies,molar_absorption,label=label)
           plot_list.append(plot)
           plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.show()
        return

    def plotAbsorption(self) :
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(self.absorption_coefficients[-1]) 
        # loop over all methods except the last, TO
        for absorptions in self.absorption_coefficients[:-1]:
            maxd = max(maxd, np.amax(absorptions) )
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        absorptions = np.array(self.absorption_coefficients[-1])
        self.absorption_coefficients[-1] = to_scale * absorptions
        plot_list = []
        plot_labels = []
        for method,shape,data,frequencies,absorptions,vf_type in zip(self.methods,self.shapes,self.shape_data,self.frequencies,self.absorption_coefficients,self.volume_fractions):
           label = method+"_"+shape+str(data)+"("+vf_type+")"
           label=label.replace(" ","")
           label=label.replace(","," ")
           label=label.replace("'","")
           label=label.replace('"','')
           plot, = pl.plot(frequencies,absorptions,label=label)
           plot_list.append(plot)
           plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.show()
        return

    def plotImaginary(self) :
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(np.imag(self.traces[-1])) 
        # loop over all methods except the last, TO
        for traces in self.traces[:-1]:
            di = np.imag(traces)
            maxd = max(maxd, np.amax(di) )
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        traces = np.array(self.traces[-1])
        self.traces[-1] = to_scale * traces
        plot_list = []
        plot_labels = []
        for method,shape,data,frequencies,traces,vf_type in zip(self.methods,self.shapes,self.shape_data,self.frequencies,self.traces,self.volume_fractions):
           label = method+"_"+shape+str(data)+"("+vf_type+")"
           label=label.replace(" ","")
           label=label.replace(","," ")
           label=label.replace("'","")
           label=label.replace('"','')
           plot, = pl.plot(frequencies,np.imag(traces),label=label)
           plot_list.append(plot)
           plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.show()
        return

    def plotReal(self) :
        """Plot loss data """
        import pylab as pl
        # Choose a suitable scale TO has to be treated differently, TO was added last
        maxd = 0.0
        max_to = 1.2*np.amax(np.real(self.traces[-1])) 
        # loop over all methods except the last, TO
        for traces in self.traces[:-1]:
            di = np.real(traces)
            maxd = max(maxd, np.amax(di) )
        # end loop over entries
        to_scale = maxd/max_to
        to_scale = 1.0
        traces = np.array(self.traces[-1])
        self.traces[-1] = to_scale * traces
        plot_list = []
        plot_labels = []
        for method,shape,data,frequencies,traces,vf_type in zip(self.methods,self.shapes,self.shape_data,self.frequencies,self.traces,self.volume_fractions):
           label = method+"_"+shape+str(data)+"("+vf_type+")"
           label=label.replace(" ","")
           label=label.replace(","," ")
           label=label.replace("'","")
           label=label.replace('"','')
           plot, = pl.plot(frequencies,np.real(traces),label=label)
           plot_list.append(plot)
           plot_labels.append(label)
        # end loop over method
        pl.legend(plot_list,  plot_labels, loc='best', numpoints=1)
        pl.show()
        return

def printReals(title, reals, no_per_line=8, format="{:9.2f}", file=sys.stdout, separator= " "):
    # 
    # Print out a list of reals prettily
    #
    len_reals = len(reals)
    print(" ",file=file)
    if title != "":
        print(title,file=file)
    nlines = int( (len_reals - 1 ) / no_per_line) + 1
    start = 0
    for i in range(nlines):
        end = start + no_per_line
        if end > len_reals: end = len_reals
        print(" "+ separator.join(format.format(r) for r in reals[start:end] ), file=file)
        start = start + no_per_line
    #end for i
    return

def print3x3(title, array, format="{:14.6f}", file=sys.stdout, separator= " "):
    # 
    # Print out a 3x3 tensor matrix
    #
    print(" ",file=file)
    if title != "":
        print(title,file=file)
    for i in range(3):
        print("      "+" ".join(format.format(p) for p in array[i]), file=file )
    # end for i
    return
