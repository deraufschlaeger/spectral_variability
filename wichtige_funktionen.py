import numpy as np
from spectres import spectres
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import pandas as pd

def list_spec_info(file_path, spec_id):

    fit=fits.open(file_path)[0]
    header = fit.header 
    imag = fit.data
    
    #Ich will erst checken wieviele Spektren wir haben
    num_of_spec = int(header['NAXIS2'])
    if not (1 <= spec_id <= num_of_spec):
            raise ValueError(f"Die Spektrum-ID muss zwischen 1 und {num_of_spec} liegen. You gave {spec_id}")

    master_str = ''
    c1 = 0
    while (header.get(f'WAT2_{c1+1:03d}')) is not None:
        master_str += header.get(f'WAT2_{c1+1:03d}') + ' '*(68-len(header.get(f'WAT2_{c1+1:03d}')))
        c1+=1

    master_list = master_str.split('spec')

    if 'multi' in master_list[0]:
        master_list.pop(0)
    if ' ' in master_list:
        master_list.remove(' ')
    if '' in master_list:
        master_list.remove('')

    spec_info = master_list[spec_id-1].split('=')[1].replace('"', '').split()
    
    if len(spec_info) <= 9:
         raise ValueError(f"es liegt nicht non-lineare Dispersion vor, was nicht implementiert ist")

    #Bennenung der Stringinformationen
    ap = float(spec_info[0])#aperture number
    beam = float(spec_info[1])#beam number
    dtype = float(spec_info[2])#dispersion type,2->non-lin, 1->log-lin, 0->lin, -1->no disp
    first_wavelength = float(spec_info[3]) 
    avg_abstand = float(spec_info[4])#lin/log->exact, non-lin->approx
    naxis = int(spec_info[5])#number of pixels
    z = float(spec_info[6])#doppler factor
    aplow = float(spec_info[7])#aperture limits
    aphigh = float(spec_info[8])
   
    w_i = float(spec_info[9])#weight
    Delta_lambda_i = float(spec_info[10])#wavelength offset
    ftype = float(spec_info[11])#type of dispersion function
    if ftype != 1:
        print(f"es liegt eine andere dispersions funktion wie chebyshev polynomial dispersion function vor -->> readmultispec")
        from readmultispec import readmultispec
        data = readmultispec(file_path)
        data1 = data['wavelen'][spec_id-1]
        data2 = data['flux'][spec_id-1]
    else:    
        order = int(spec_info[12]) #Ordnungen der Chebyshev Polynome
        p_min = float(spec_info[13])
        p_max = float(spec_info[14])
        coeff = np.asarray(spec_info[15:], dtype=float) # Coefficients der Ordnungen
        if len(coeff) != order:
            raise AttributeError(f'Error while Reading Data: number of coefficients doesnt fit the order number')
    
        pmiddle = (p_max+p_min)/2
        prange = p_max-p_min
        n = (np.arange(naxis, dtype=float) + 1 - pmiddle)/(prange/2) 
    
        x0 = np.ones(naxis, dtype=float)
        x1 = n
        wave = x0 * coeff[0] + x1 * coeff[1]
        for i in range(2, order):
            x2 = 2 * n * x1 - x0
            wave += x2 * coeff[i]
            x0 = x1
            x1 = x2
    
        wavelengths = w_i*(Delta_lambda_i+wave)/(1+z)
        
        pixel_val = imag[spec_id-1,]
    
        data1 = [wavelengths[i] for i in range(len(wavelengths))]
        data2 = [pixel_val[i] for i in range(len(wavelengths))]
    
        
        #for index, i in enumerate(data2):
            #if i < 0 :
                #data2[index] = 0
                
    return data1, data2

def get_master_lists(file_path):
    fit=fits.open(file_path)[0]
    header = fit.header 
    
    #Ich will erst checken wieviele Spektren wir haben
    num_of_spec= int(header['NAXIS2'])
    master_wave= []
    master_flux = []
    for i in range(1, num_of_spec+1):
        wave, flux = list_spec_info(file_path, i)
        master_wave.append(wave)
        master_flux.append(flux)
    return master_wave, master_flux

def calculate_bin_averages(bins, x_values, y_values):
    # Ergebnisliste für die Durchschnittswerte der Bins
    average_x_values = []
    average_y_values = []

    for bin in bins:
        start_index, end_index = bin

        if end_index-start_index <= 50:
            # Hole die entsprechenden Werte aus den Arrays
            x_bin_values = x_values[start_index:end_index+1]
            y_bin_values = y_values[start_index:end_index+1]
    
            # Berechne den Durchschnitt für die X- und Y-Werte
            avg_x = sum(x_bin_values) / len(x_bin_values)
            avg_y = sum(y_bin_values) / len(y_bin_values)
            
            # Füge die Durchschnittswerte zur Ergebnisliste hinzu
            average_x_values.append(avg_x)
            average_y_values.append(avg_y)

        else:
            # teile den bin in zwei teile
            middle = int((start_index+end_index)/2)
            
            # Hole die entsprechenden Werte aus den Arrays
            x_bin_values_1 = x_values[start_index:middle+1]
            x_bin_values_2 = x_values[middle:end_index+1]
            y_bin_values_1 = y_values[start_index:middle+1]
            y_bin_values_2 = y_values[middle:end_index+1]
            
            # Berechne den Durchschnitt für die X- und Y-Werte
            avg_x_1 = sum(x_bin_values_1) / len(x_bin_values_1)
            avg_x_2 = sum(x_bin_values_2) / len(x_bin_values_2)
            avg_y_1 = sum(y_bin_values_1) / len(y_bin_values_1)
            avg_y_2 = sum(y_bin_values_2) / len(y_bin_values_2)
            
            # Füge die Durchschnittswerte zur Ergebnisliste hinzu
            average_x_values.append(avg_x_1)
            average_x_values.append(avg_x_2)
            average_y_values.append(avg_y_1)
            average_y_values.append(avg_y_2)
    
    return average_x_values, average_y_values

def Blaze_estimate_bins(x, y, spec_index):
    spec_id = spec_index + 1
    ''' bins_x = []
    bins_y = []'''
    
    
    if spec_id == 1:#rauschen,sonst ok
        bin_indices = [[3, 20], [50, 100], [150, 170], [220, 250], [300, 400], [550, 600], [695, 725], [840, 870], [900, 950], [1220, 1270], [1310, 1380], [1450, 1500], [1550, 1600], [1660, 1800], [1850, 1900], [2000, 2030]]
    
    if spec_id == 2:#!
        bin_indices = [[3, 20], [30, 50], [100, 150], [175, 225], [300, 350], [425, 475], [650, 851], [1100, 1150], [1350, 1400], [1650, 1700], [1780, 1800], [1950, 2000], [2020, 2040]]
    
    if spec_id == 3:#!
        bin_indices = [[5, 20], [50, 100], [150, 200], [320, 360], [600, 650], [820, 860], [1150, 1200], [1310, 1360], [1550, 1580], [1750, 1800], [1900, 1950], [2020, 2040]]
    
    if spec_id == 4:#ok
        bin_indices = [[2, 20], [70, 120], [500, 550], [750, 800], [950, 1000], [1100, 1150], [1280, 1330], [1600, 1700], [1830, 1850], [1980, 2000], [2030, 2045]]
    
    if spec_id == 5:#ok
        bin_indices = [[5, 70], [90, 140], [200, 250], [400, 450], [550, 600], [865, 910], [1060, 1110], [1300, 1350], [1500, 1550], [1670, 1710], [1850, 1900], [2030, 2045]]
    
    if spec_id == 6:#ok
        bin_indices =[[3, 20], [80, 130], [300, 350], [580, 630], [750, 800], [1030, 1080], [1250, 1300], [1380, 1390], [1530, 1580], [1610, 1650], [1850, 1900], [1975, 2000], [2020, 2040]]
    
    if spec_id == 7:#cheby nicht gut
        bin_indices = [[3,13], [50, 70], [310, 340], [550, 600], [1220, 1270], [1400, 1450], [1650, 1700], [1750, 1800], [1920, 1970], [2020, 2040]]
    
    if spec_id == 8:#ok
        bin_indices = [[5, 20], [150, 200], [220, 270], [500, 550], [700, 750], [800, 820], [930, 980], [1180, 1230], [1560, 1610], [1750, 1790], [1950, 2000], [2030, 2040]]
    
    if spec_id == 9:#ok
        bin_indices = [[5, 15], [30, 80], [140, 180], [550, 600], [700, 750], [830, 860], [1155, 1180], [1270, 1320], [1450, 1500], [1530, 1560], [1580, 1610], [1750, 1790], [1900, 1950], [2043, 2047]]
    
    if spec_id == 10:#ok
        bin_indices = [[3, 17], [70, 90], [200, 250], [440, 470], [600, 650], [760, 810], [950, 1000], [1200, 1250], [1420, 1460], [1705, 1750], [1950, 2000], [2020,2040]]
    
    if spec_id == 11:#ok
        bin_indices = [[3, 7], [10, 60], [200, 250], [310, 350], [550, 600], [735, 755], [850, 900], [1230, 1280], [1370, 1410], [1550, 1600], [1700, 1750], [1850, 1890], [1950, 2000], [2030, 2042]]
    
    if spec_id == 12:#ok
        bin_indices = [[4, 14], [20,200], [270, 305], [540, 580], [830, 850], [1000, 1050], [1180, 1230], [1290, 1330], [1450, 1500], [1750, 1800], [1900, 1950], [2010, 2025], [2035, 2045]]
    
    if spec_id == 13:#cheby nicht gut
        bin_indices = [[5, 15], [50, 65], [150, 200], [250, 300], [440, 460], [640, 690], [900, 950], [1190, 1220], [1420, 1440], [1550, 1570], [1850, 1900], [2030, 2040]]
    
    if spec_id == 14:#mit peak 4100 im vega spektrum vergleichen
        bin_indices = [[5, 15], [60, 90], [280, 320], [460, 500], [650,750], [950, 1000], [1160, 1190], [1350, 1400], [1550, 1580], [1650, 1700], [1770, 1800], [1810, 1860], [1970, 2020], [2040, 2045]]
    
    if spec_id == 15:
        bin_indices = [[4, 14], [50, 70], [120, 170], [320, 370], [450, 500], [720, 970], [1070, 1120], [1355, 1500], [1490, 1540], [1600, 1650], [1750, 1780], [1850, 1900], [1980, 2010], [2030, 2045]]
    
    if spec_id == 16:
        bin_indices = [[3, 15], [30, 80], [150, 200], [380, 400], [480, 520], [650, 700], [890, 940], [1300, 1350], [1420, 1470], [1550, 1650], [1780, 1900], [1965, 2000], [2040, 2047]]
    
    if spec_id == 17:#ok
        bin_indices = [[3,13], [50, 80], [190, 220], [300, 320], [380, 400], [450, 480], [635, 663], [750, 800], [900, 950], [1090, 1130], [1250, 1300], [1550, 1570], [1650, 1700], [1760, 1800], [1900, 1950], [2030, 2045]]
    
    if spec_id == 18:#ok~
        bin_indices = [[5, 10], [100, 150], [250, 300], [550, 600], [755, 790], [1000, 1050], [1250, 1300], [1500, 1600], [1725, 1750], [1780, 1800], [1950, 2000], [2020, 2047]]
    
    if spec_id == 19:#ok
        bin_indices = [[0, 3], [30, 50], [150, 190], [220, 280], [420, 450], [535, 650], [700, 750], [1035, 1080], [1250, 1290], [1464, 1503], [1640, 1670], [1800, 1950], [2020, 2030]]
    
    if spec_id == 20:#ok vega spektrum vergleichen
        bin_indices = [[3, 40], [50, 100], [250, 300], [350, 380], [540, 575],  [650, 700], [850, 900], [1050, 1100], [1300, 1350], [1450, 1500], [1600, 1650], [1700, 1750], [1890, 1907], [2030, 2045]]
    
    if spec_id == 21:#cheby nicht gut, sonst ok
        bin_indices = [[5, 30], [50, 90], [200, 300], [400, 440], [600,650], [800, 850], [970, 1020], [1300, 1330], [1575, 1625], [1750, 1800], [1920, 1970], [2040, 2045]]
        bin_indices_2 = [[5, 40], [50, 90], [200, 300], [350, 440], [500, 550], [1400, 1450], [1595, 1660], [1780, 1800], [1810, 1880], [1920, 1970], [1990, 2030]]
    
    if spec_id == 22:#ok
        bin_indices = [[3, 13], [70, 100], [170, 180], [250, 300], [400, 550], [640, 680], [840, 890], [1180, 1230], [1400, 1450], [1580, 1620], [1810, 1860], [1900, 1930], [2030, 2044]]
    
    if spec_id == 23:#ok
        bin_indices = [[3, 13], [120, 170], [300, 350], [550, 600], [730, 760], [1050, 1100], [1200, 1250], [1400, 1450], [1500, 1550], [1700, 1750], [1880, 1900], [1960, 2020], [2030, 2044]]
    
    if spec_id == 24:#ok
        bin_indices = [[0, 5], [50, 125], [300, 350], [500, 550], [680, 730], [865, 885], [1050, 1100], [1230, 1270], [1420, 1460], [1580, 1620], [1700, 1730], [1900, 1920], [2042, 2045]]
    
    if spec_id == 25:#ok
        bin_indices = [[3, 20], [50, 75], [150, 200], [310, 355], [560, 580], [700, 750], [880, 930], [1080, 1127], [1270, 1300], [1460, 1510], [1600, 1650], [1800, 1850], [2000, 2018], [2020, 2045]]
    
    if spec_id == 26:#ok
        bin_indices = [[2, 23], [63, 110], [270, 300], [340, 370], [480, 510], [600, 750], [883, 917], [1065, 1100], [1300, 1350], [1520, 1565], [1650, 1690], [1850, 1900], [2010, 2045]]
    
    if spec_id == 27:#ok
        bin_indices = [[3, 13], [20, 100], [150, 173], [210, 230], [420, 470], [513, 675], [740, 780], [940, 970], [1115, 1135], [1385, 1405], [1590, 1610], [1777, 1810], [1990, 2020], [2035, 2046]]
    
    if spec_id == 28:#ok
        bin_indices = [[3, 13], [70, 90], [190, 210], [280, 320], [510, 540], [705, 715], [865, 875], [995, 1004], [1180, 1200], [1325, 1345], [1565, 1575], [1720, 1770], [1980, 2010], [2040, 2045]]
    
    if spec_id == 29:#ok~
        bin_indices = [[3, 20], [130, 150], [200, 250], [400, 420], [640, 680], [800, 820], [990, 1010], [1300, 1350], [1550, 1600], [1700, 1750], [1880, 1950], [2000, 2020], [2030, 2045]]
    
    if spec_id == 30:#ok
        bin_indices = [[3, 7], [90, 177], [240, 290], [350, 400], [500, 550], [720, 770], [950, 1000], [1100, 1150], [1285, 1330], [1420, 1460], [1550, 1600], [1750, 1800], [1900, 1950], [2030, 2046]]
    
    if spec_id == 31:#ok
        bin_indices = [[3, 10], [30, 50], [125, 150], [200, 250], [310, 330], [500, 550], [650, 970], [1020, 1070], [1320, 1340], [1460, 1500], [1590, 1640], [1705, 1750], [1800, 1900], [1950, 2000], [2030, 2045]]
    
    if spec_id == 32:#ok, cheby ~
        bin_indices = [[3, 20], [100, 150], [360, 380], [500, 550], [700, 750], [950, 1000], [1150, 1180], [1400, 1500], [1560, 1610], [1740, 1760], [1810, 2000], [2020, 2040]]
    
    if spec_id == 33:#ok
        bin_indices = [[3, 20], [30, 50], [240, 280], [400, 450], [550, 600], [740, 785], [1010, 1060], [1100, 1250], [1410, 1460], [1880, 1930], [2030, 2047]]
        
    if spec_id == 34:#ok
        bin_indices = [[3, 30], [50, 100], [380, 400], [785, 810], [1050, 1100], [1300, 1330], [1395, 1445], [1630, 1650], [1700, 1950], [2040, 2047]]
    
    if spec_id == 35: #n@j@!
        bin_indices = [[30, 50], [120, 170], [230, 250], [350, 400], [660, 710], [800, 850], [1055, 1065], [1300, 1350], [1600, 1650], [1750, 1950], [2000, 2045]]
    # [640, 700], [1130, 1190],
    if spec_id == 36:#ok
        bin_indices = [[3, 10], [40, 60], [220, 260], [370, 550], [700, 750], [1080, 1100], [1300, 1320], [1500, 1510], [1630, 1680], [1710, 2000], [2020, 2040]]
    
    if spec_id == 37:#mh
        bin_indices = [[3, 10], [85, 95], [220, 250], [400, 450], [615, 640], [785, 805], [930, 970], [1115, 1130], [1727, 1743], [1920, 1930], [2023, 2027]]
    
    if spec_id == 38:
        bin_indices = [[3,15], [30, 50], [170, 270], [420, 470], [550, 670], [750, 800], [850, 900], [1130, 1170], [1330, 1370], [1450, 1490], [1710, 1760], [1850, 1890], [1970, 1980], [2040, 2045]]
    
    if spec_id == 39: #mh
        bin_indices = [[3, 20], [40, 90], [190, 230], [248, 351], [500, 530], [800, 850], [1000, 1040], [1150, 1190], [1300, 1350], [1450, 1490], [1690, 1740], [1900, 1950], [1970, 2000], [2035, 2037]]
    
    if spec_id == 40:#ok
        bin_indices = [[3, 20], [100, 140], [320, 370], [430, 480], [550, 600], [800, 850], [1000, 1050], [1095, 1140], [1210, 1250], [1360, 1400], [1555, 1600], [1720, 1740], [1905, 1950], [2035, 2046]]
    
    if spec_id == 41:# am rechtem rand komisch aber ok
        bin_indices = [[15, 20], [70, 120], [170, 220], [250, 300], [390, 420], [450, 500], [650, 700], [985, 1005], [1460, 1490], [1710, 1740], [1860, 1890], [1950, 1970], [2020, 2030], [2043, 2045]]
    
    if spec_id == 42:
        bin_indices = [[3, 5], [140, 155], [250, 290], [390, 440], [480, 510], [560, 610], [640, 680], [1090, 1120], [1480, 1530], [1650, 1660], [1810, 1840], [1975, 1990], [2030, 2043]]
    
    if spec_id == 43:#ok
        bin_indices = [[2, 10], [80, 100], [260, 275], [470, 500], [550, 570], [760, 800], [1150, 1200], [1420, 1460], [1510, 1520], [1900, 1960], [2030, 2040]]
    
    if spec_id == 44:#ok
        bin_indices = [[3, 10], [60, 100], [210, 240], [450, 480],[680, 730], [830, 860], [1050, 1070], [1545, 1565], [1750, 1780], [1940, 1955], [2030, 2045]]
    
    if spec_id == 45:#ok
        bin_indices = [[3, 8], [170, 220], [600, 620], [760, 780], [1285, 1330], [1460, 1490], [1600, 1620], [1800, 1815], [1930, 1945], [2010, 2025]]
    
    if spec_id == 46:#ok
        bin_indices = [[3, 10], [95, 105], [190, 200], [300, 330], [505, 525], [580, 600], [800, 850], [1190, 1220], [1520, 1540], [1750, 1800], [1885, 1905], [2000, 2020], [2040, 2045]]
    
    if spec_id == 47:
        bin_indices = [[5, 15], [20, 70], [90, 110], [120, 150], [185, 235], [290, 450], [550, 600], [700, 750], [820, 880], [920, 970], [1045, 1150], [1170, 1300], [1370, 1460], [1520, 1560], [1570, 1630], [1640, 1710], [1740, 1810], [1870, 1960], [2030, 2045]]
        
    x_spline, y_spline = calculate_bin_averages(bin_indices, x, y)

    #print(x_spline)

    cs = CubicSpline(np.array(x_spline), np.array(y_spline))
    return cs
    
def delta(flux, neighbor):
    flux = list(flux)
    lower = [flux[0]] * neighbor
    upper = [flux[-1]] * neighbor
    master = lower + flux + upper
    out = []
    for f in flux:
        out.append(f-0.5*(master[master.index(f, neighbor, -neighbor)-neighbor] + master[master.index(f, neighbor, -neighbor)+neighbor]))
    out = np.array(out)
    return out

def gauss_poly(x, A, mu, sigma, a0):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) + a0# + a1 * x + a2 * x ** 2

def remove_outliers(array, x0, sigma):
    out = []
    for i in array:
        if x0-2*sigma <= i <= x0+2*sigma:
            out.append(i)
    return out

def get_sorted_list_and_per(list1, per=True):
    L = [(list1[i],i) for i in range(len(list1))]
    L.sort()
    sorted_list, permutation = zip(*L)
    sorted_list = list(sorted_list)
    if per == True:
        return sorted_list, permutation
    if per == False:
        return sorted_list

def sort_with_per(to_sort_list, permutation):
    sorted_list = []
    for index2 in permutation:
        sorted_list.append(to_sort_list[index2])
    return sorted_list

def get_signal_to_noise(flux, **kwargs):
    neighbor = kwargs.get('neighbor', 2)
    conf_level = kwargs.get('conf_level', 0)
    if conf_level != -1 and conf_level != 0 and conf_level != 1 and conf_level != 2:
        raise ValueError("Usage error in get_signal_to_noise: Invalid value for qualifier 'conf_level'")
    ndist = delta(flux, neighbor)
    if conf_level != -1:
        key_list = ['noise_min', 'noise', 'noise_max', 'snr_min', 'snr', 'snr_max']
    else:
        key_list = ['noise', 'snr']
    out = dict.fromkeys(key_list, [])
    if len(np.nonzero(ndist)) == 0:
        out['noise'] = 0
        out['snr'] = np.inf
        if conf_level != -1:
            out['noise_min'] = 0
            out['noise_max'] = 0
            out['snr_min'] = np.inf
            out['snr_max'] = np.inf
    else:
        ndist_filtered = remove_outliers(ndist, np.average(ndist), np.std(ndist))
        hist, bin_edges = np.histogram(ndist_filtered, bins=200)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist_errors = np.sqrt(hist)
        hist_errors[hist_errors == 0] = 1

        #pos_shift = abs(2 * bin_centers[0])
        #bin_centers += pos_shift

        #obs = define_couts(hist)

        p0 = [len(ndist_filtered), np.average(ndist_filtered), np.std(ndist_filtered), 0]
        bounds = ([-np.inf, -np.inf, 0.1*np.std(ndist_filtered), -np.inf], [np.inf, np.inf, 10*np.std(ndist_filtered)+10**-9, np.inf])

        para, cov = curve_fit(gauss_poly, bin_centers, hist, p0=p0, bounds = bounds, sigma = hist_errors)
#, sigma = hist_errors
        conv_factor = np.sqrt(1.5)
        f_med  = np.median(flux)
        noise = para[2]/conv_factor
        out['noise'] = noise
        out['snr'] = f_med/noise
        if conf_level != -1:
            perr = np.sqrt(np.diag(cov))
            noise_min = (para[2]-perr[2])/conv_factor
            out['noise_min'] = noise_min
            noise_max = (para[2]+perr[2])/conv_factor
            out['noise_max'] = noise_max
            out['snr_min'] = f_med / noise_min
            out['snr_max'] = f_med / noise_max
            
    return out

def optimize_wavegrid_echelle(l, f, res):
    lambda_optimal=[]
    for i in range(len(l)):
        lambda_optimal.append(l[0]*((4*res+1)/(4*res-1))**i)
    sorted_l, per = get_sorted_list_and_per(l)
    sorted_f = sort_with_per(f, per)
    
    flux_optimal = spectres(np.array(lambda_optimal), np.array(sorted_l), np.array(sorted_f), fill = 0, verbose=False)
    lambda_optimal = np.array(lambda_optimal)
    return lambda_optimal, flux_optimal

def get_signal_to_noise_curve(l, f, data_points=3000, conf_level=0):
    len_l = len(l)
    if conf_level != -1:
        out = {'l':[], 'noise_min':[],'noise':[],'noise_max':[],'snr_min':[],'snr':[],'snr_max':[]}
    else:
        out = {'l':[],'noise':[],'snr':[]}
    index = 0
    if data_points > len_l:
        data_points = len_l
    data_points -= 1 #conversion to indices
    while index + data_points <= len_l - 1:
        #slice = [index:index+data_points:1]
        flux_temp = f[index:index+data_points]
        out_temp = get_signal_to_noise(flux_temp, conf_level = conf_level)
        for key in out_temp.keys():
            out[key].append(out_temp.get(key))
        out['l'].append(l[int(index+data_points/2)])
        index += int(data_points/2)
    #the last interval needs special treatment
    #slice =  [len-1-data_points : len-1 : 1]
    flux_temp = f[len_l-1-data_points:len_l-1]
    out_temp = get_signal_to_noise(flux_temp, conf_level = conf_level)
    for key in out_temp.keys():
        if key not in out or not isinstance(out[key], list):
            out[key] = []
        out[key].append(out_temp.get(key))
        out[key] = np.array(out[key])
    out['l'].append(l[int(len_l-data_points/2)])
    out['l'] = np.array(out['l'])
    return out

def coadd(file_path):
    header = fits.open(file_path)[0].header
    n_orders = header['NAXIS2']
    slit = header['CO-SLIT']
    l_in = [list_spec_info(file_path, i)[0] for i in range(1,n_orders+1)]
    f_in = [list_spec_info(file_path, i)[1] for i in range(1,n_orders+1)]
    u = []
    for index, order_l in enumerate(l_in):
        u.extend(order_l)
        f_in[index] /= np.median(f_in[index])
    u = np.array(sorted(u))
        
    if slit.strip() == '0.90mm':
        res = 35000 #source: kso.tls-tautenburg.de/TLS/index.php?id=31 // nachgeprüft
    elif slit.strip() == '1.00mm' or slit.strip() == '1.0mm':
        res = 32000 # nachgeprüft
    elif slit.strip() == '0.52mm':
        res = 50000
    else:
         raise ValueError(f"Spaltbreite von {slit} muss noch implementiert werden")
        
    l_new = [optimize_wavegrid_echelle(u, np.interp(u, l_in[i],f_in[i], left = 0, right = 0), res)[0] for i in range(n_orders)]
    f_new = [optimize_wavegrid_echelle(u, np.interp(u, l_in[i],f_in[i], left = 0, right = 0), res)[1] for i in range(n_orders)]

    len2 = len(l_new[0])
    f = []
    norm = []
    for i in range(n_orders):
        l_new[i] = l_new[i][f_new[i]!=0]
        f_new[i] = f_new[i][f_new[i]!=0]
        
        snr = get_signal_to_noise_curve(l_new[i], f_new[i], data_points = 200, conf_level = -1 )
        weight = (np.interp(l_new[i], snr['l'], snr['snr']))**2
        norm.append(weight)
        f.append(weight*f_new[i])
    for index, li in enumerate(norm):
        f[index] /= li

    return l_new, f

def coadd_list(l_in, f_in, file_path):
    header = fits.open(file_path)[0].header
    n_orders = header['NAXIS2']
    slit = header['CO-SLIT']
    
    u = []
    for index, order_l in enumerate(l_in):
        u.extend(order_l)
        f_in[index] /= np.median(f_in[index])
    u = np.array(sorted(u))
        
    if slit.strip() == '0.90mm':
        res = 35000 #source: kso.tls-tautenburg.de/TLS/index.php?id=31 // nachgeprüft
    elif slit.strip() == '1.00mm' or slit.strip() == '1.0mm':
        res = 32000 # nachgeprüft
    elif slit.strip() == '0.52mm':
        res = 50000
    else:
         raise ValueError(f"Spaltbreite von {slit} muss noch implementiert werden")
        
    l_new = [optimize_wavegrid_echelle(u, np.interp(u, l_in[i],f_in[i], left = 0, right = 0), res)[0] for i in range(n_orders)]
    f_new = [optimize_wavegrid_echelle(u, np.interp(u, l_in[i],f_in[i], left = 0, right = 0), res)[1] for i in range(n_orders)]

    len2 = len(l_new[0])
    f = []
    norm = []
    for i in range(n_orders):
        l_new[i] = l_new[i][f_new[i]!=0]
        f_new[i] = f_new[i][f_new[i]!=0]
        
        snr = get_signal_to_noise_curve(l_new[i], f_new[i], data_points = 200, conf_level = -1 )
        weight = (np.interp(l_new[i], snr['l'], snr['snr']))**2
        norm.append(weight)
        f.append(weight*f_new[i])
    for index, li in enumerate(norm):
        f[index] /= li

    return l_new, f

def group_average(group):
    flux = np.array(group["flux"])
    snr = np.where(np.array(group["snr"]) == 0,1,np.array(group["snr"])*gauss_poly(group["index_in_order"], 1, group["len_of_order"]/2, group["len_of_order"]/10, 0))
    return np.average(flux, weights = snr)

def weighted_average(group):
    # Berechne den gewichteten flux
    weighted_flux = group_average(group)
    
    # Optional: Berechnung des gewichteten snr (oder andere Statistiken)
    avg_snr = group['snr'].mean()
    
    return pd.Series({'weighted_flux': weighted_flux, 'avg_snr': avg_snr})

def calc_snr_pointwise(intensity, bin_size):
    snr = np.ones(len(intensity))
    for i in range(len(intensity)):
        bin_slice = intensity[i:i + bin_size]
        if np.std(bin_slice) != 0:
            snr[i] = np.mean(bin_slice) / np.std(bin_slice)
        else:
            snr[i] = 0
    return snr

def merge_orders(file_path):
    """
    Kombiniert mehrere Spektren aus einer FITS-Datei, um einen gewichteten Durchschnitt von 
    Flux und Wellenlängen zu berechnen.

    Diese Funktion liest Spektren aus der angegebenen FITS-Datei und 
    berechnet für jeden Wellenlängenwert den gewichteten Durchschnitt des Flux. 
    Zusätzlich werden Signal zu Rausch Verhältnisse (SNR) für die Fluxwerte berechnet. Vor der 
    Berechnung wird der Flux durch eine Schätzung des Blaze-Profils normalisiert.

    Parameters
    ----------
    file_path : str
        Der Pfad zur fertig reduzierten FITS-Datei, die die Spektrendaten enthält.

    Returns
    -------
    tuple
        Ein Tupel, das die Wellenlängen und die entsprechenden gewichteten Fluxwerte enthält.
        - wavelengths : list
            Eine Liste der Wellenlängenwerte.
        - weighted_flux : list
            Eine Liste der gewichteten Fluxwerte, die für jede Wellenlänge berechnet wurden.

    Notes
    -----
    Die Funktion verwendet eine Bin Größe von 20 für die SNR-Berechnung. 
    Ungültige Fluxwerte (z.B. gleich 0) werden entfernt, und die Wellenlängen- und Fluxarrays 
    werden um die Randwerte gekürzt, um die Datenqualität zu verbessern.

    Example
    --------
    >>> wavelengths, fluxes = merge_orders('path/to/spectrum_file.fits')
    """
    bin_size = 20

    #list of lists of spectral orders
    master_wave, master_flux = get_master_lists(file_path)
    
    #BLAZE-FUNKTION
    for i in range(len(master_flux)):
        master_flux[i] /= Blaze_estimate_bins(master_wave[i], master_flux[i], i)(master_wave[i])
    #BLAZE-FUNKTION 
    
    master_wave, master_flux = coadd_list(master_wave, master_flux, file_path)
    
    for i in range(len(master_wave)):
        master_wave[i] = master_wave[i][master_flux[i] != 0]
        master_flux[i] = master_flux[i][master_flux[i] != 0]
        master_wave[i] = master_wave[i][1:-1]
        master_flux[i] = master_flux[i][1:-1]

    data = {
        "wavelength": [],
        "flux": [],
        "snr": [],
        "order_index": [],
        "index_in_order": [],
        "len_of_order": []
    }
    
    # Sammle alle Wellenlängen, fluxe und SNR in einem DataFrame
    for order_index, (order_wave, order_flux) in enumerate(zip(master_wave, master_flux)):
        snr = calc_snr_pointwise(order_flux, bin_size)
        data["wavelength"].extend(order_wave)
        data["flux"].extend(order_flux)
        data["snr"].extend(snr)
        data["order_index"].extend([order_index] * len(order_wave))
        data["index_in_order"].extend(list(range(len(order_wave))))
        data["len_of_order"].extend([len(order_wave)] * len(order_wave))
        
    df = pd.DataFrame(data)
    combined_flux = []
    combined_wave = []
        
        # Gruppiere nach Wellenlänge, um den flux von mehrfach vorkommenden Wellenlängen zu berechnen
    grouped = df.groupby(["wavelength"], as_index=False).apply(weighted_average, include_groups=False).reset_index()
    return grouped["wavelength"], grouped["weighted_flux"]
