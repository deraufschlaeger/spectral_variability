procedure red_echelle_VIS (inimages, outimages, tharimages)

#   Das Programm macht die ersten Datenreduktions-
#   schritte fuer die Echelle-Daten aus Tautenburg.
#
#   Dies ist die neue Version fuer den ganz neuen 2kx2k Andor 
#   Chip, der ab 2021 verwendet wird. Die alte Version war
#   fuer den Chip der von 2015 bos 2021 verwendet wurde
#   
#   Es ist identisch zu redech_frame, macht aber
#   ausserdem noch die Wellenlaengenkalibration.
#
#neue Wellenlaengenkalibration: apall vega.fits (alles yes), apall Thar.fits (nur extract: yes)
#   Folgende Frames werden benoetigt:
#
#   Slope.fits   : Struktur des Darks, Bias abgezogen.
#   [2048,2048]   Die Werte sind also alle etwa Null.
#
#   Flat.fits    : Flatfield image, das in den Ordnungen
#   [2048,2048]   nur Werte um 1.0 enthaelt. Das
#                 Wird einfach dadurch hergestellt,
#                 dass man alle Flatfields nimmt, mittelt,
#                 den Bias abzieht und dann das Zeug durch
#                 apflatten jagt.
#
#   find_orders : Referenzframe fuer die Ordnungen.
#   [2048,2048]   Das ist einfach die Summe aller
#                 SCI-Frames. Dann einmal durch apall
#                 geschickt, um Parameter festzulegen.
#
#   ThAr        : ThAr-image, in dem die Linien
#   [2048,51]     bereits identifiziert wurden.
#    
#   ------------------------------------------------------------
#
#   Vorbereitung fuer VIS Channel:
#
#   imcopy Flat_vis.fits Flat.fits 
#   imcopy ThAr_vis.fits ThAr.fits
#   imcopy find_orders_vis.fits find_orders.fits
#   apall find_orders.fits
#   dispcor : table  = wavetable_ir
#   ecreidentify: referenc = ThArXI_ir,
#     oder Frame ThAr einmal durch ecreidentify laufen lassen.
#   echelle
#      dispaxis=1; verbose=no
#   epar refspec, apall, apscatter, ecreidentify (thresho 300)
#
#   ------------------------------------------------------------
#
#   Hinweise:
#      
#   Das Format der CCDs hate sich im Mai geandert.
#   jetzt ist es nur noch 2147,2063 statt 2148,2063
#           
#   Um das richtige Format zu erhalten:
#   imcopy frame.fits[1:2048,1:2048] frame1.fits
#   Die Ordnungen muessen immer von rechts nach links benannt
#   werden, die groesste Wellenlaenge erhaelt die Ordnung 1.
#
#   Die Frames brauchen nur ent-cosmic-ed zu werden 
#
#   inimages   : raw image
#   outimages  : output image
#   tharimages : das dazugehoerige thar
#
#   redech_thar in.fits out.fits thar.fits
#
#  17.7.2021                                         Eike Guenther

string inimages   {prompt="input image"}
string outimages  {prompt="output image"}
string tharimages {prompt="thar (raw-frame)"}
struct *list1
struct *list2
struct *list3

begin
     file infile
     file outfile
     file tharfile
     string in
     string out
     string thar
 
# --- Erzeugung von temporaeren Filenamen:
      infile   = mktemp ("tmp")
      outfile  = mktemp ("tmp")
      tharfile = mktemp ("tmp")
 
# --- Umwandeln der Listen von Frames in temporaere Files:
      sections(inimage,   option="root", > infile)
      sections(outimage,  option="root", > outfile)
      sections(tharimage, option="root", > tharfile)
      list1 = infile
      list2 = outfile
      list3 = tharfile
 
# --- Bearbeiten der einzelnen Frames in einer Schleife
 
      print("********************************************")
      while (fscan(list1,in) !=EOF){
 
          if (fscan(list2, out) == EOF){
              print (" Not enough output frames ")
              return
              }
          if (fscan(list3, thar) == EOF){
              print (" Not enough thar frames ")
              return
              }

# --- Loeschen der dummies :
      print("********************************************")
      print(" ")
      print("     Automatic data reduction of the")
      print("          TLS Echelle spectra")
      print("       taken with the Andor camera")
      print(" ")
      print("             Version 3.0")
      print(" ")
      print("********************************************")
      imdelete(out, verify=yes)
      print("     Current Frame:     ")
      imhead(in,longheader=no)
      print("********************************************")

# --- Start der Berechungen :
      imcopy(in,"dumaa")
#      print ("removing cosmic rays [this will take some time]")
#      cosmicrays(input="dumaa",output="duma",thresho=100,fluxrat=0.5,npasses=14,window=7,interac=no)
      imcopy("dumaa","duma")
      print ("cosmic rays not removed")
      imdelete("dumaa",verify=no)
      hedit("duma",fields="DISPAXIS",value=1,add=yes)
      hedit("duma",fields="INSTRUMENT",value="TLS-Echelle",add=yes)
      print("removing slope")
      imarith("duma","-","Slope","dumc",pixtype="real",calctyp="real")
      imdelete("duma",verify=no)
      print ("removing the BIAS offset:")
      imarith("dumc","-",300.9,"dumd",pixtype="real",calctyp="real")
#      ccdproc("dumc",output="dumd",fixpix=no, overscan=yes, biassec="[1:2,2:2046]")
      imdelete("dumc",verify=no)
      print("extracting image")
      imcopy("dumd","dume") 
      imdelete("dumd",verify=no)
      print("flat-fielding the image")
      imarith("dume","/","flat_VIS","dumf"+"_"+in,pixtype="real",calctyp="real")      
      imdelete("dume",verify=no)
      print(" ")
      print("removing scattered light")
      print(" ")
      apscatter("dumf"+"_"+in,"dumi",referen="find_orders", line=1024, find=yes, recente=yes)
      print (" ")
      print("extracting the orders")
      print(" ")
      apall("dumi",output="dumj",references="dumf"+"_"+in,backgro="none",trace=no, resize=no)
      imdelete("dumi",verify=no)
      print ("")
      print("extracting the thar")
      imcopy(thar,"dumta")
      hedit("dumta",fields="DISPAXIS",value=1,add=yes)
      hedit ("dumta",fields="DATASEC", add=no, delete=yes)
      hedit ("dumta",fields="CCDSEC", add=no, delete=yes)
      hedit ("dumta",fields="BIASSEC", add=no, delete=yes)
      hedit ("dumta",fields="BIASSEC1", add=no, delete=yes)
      imarith("dumta","-",300.9,"dumtb",pixtype="real",calctyp="real")
#      ccdproc("dumta",output="dumtb",fixpix=no, overscan=yes, biassec="[1:2,2:2046]")
      imdelete("dumta",verify=no)
      imcopy("dumtb","dumtc")
      imdelete("dumtb",verify=no)
      print(" ")
      apall("dumtc",output="dumtd",interac=no,references="dumf"+"_"+in,backgro="none")
      imdelete("dumf"+"_"+in,verify=no)  
      imdelete("dumtc",verify=no)
      ecreidentify("dumtd","ThAr")
      refspectra("dumj",references="dumtd",sort="",group="")
      dispcor("dumj","dumk",linearize=no)
      imdelete("dumtd",verify=no)
    
#     Trick: Die Werte des nicht-rebinnten Spektrums dumj werden       
#     in das File mit dem rebinnten Spektrum (dumk) eingetragen.
      imreplace("dumk",0.0,lower=INDEF,upper=INDEF)
      imarith ("dumk","+","dumj",out,pixtype="real",calctyp="real")
      imdelete("dumj",verify=no)
      imdelete("dumk",verify=no)
      print("********************************************")
      print("Wavelength calibrated spectrum created ")
      imhead(out,longheader=no)
      print("********************************************")
      }

# --- Saubermachen
      list1= ""
      list2= ""
      list3= ""
      delete (infile, verify=no)
      delete (outfile, verify=no)
      delete (tharfile, verify=no)
end
