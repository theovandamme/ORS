import ee
import colorsys
import numpy as np

def change(pre, post, index): 
    # Compute the specified index for both pre and post images
    index_pre = pre.map(globals()[index])
    index_post = post.map(globals()[index])
    
    if index in ['HSV_1', 'HSV']:
        # Create mosaics for the HSV indices
        index_pre = index_pre.qualityMosaic(index)
        index_post = index_post.qualityMosaic(index)
        # Calculate relative difference index for HSV
        rd_index = ((index_post.subtract(index_pre)).divide((index_post.add(index_pre)))).multiply(100)
    else: 
        # Create quality mosaics for other indices
        index_pre = index_pre.qualityMosaic(index).rename(index + '_pre')
        index_post = index_post.qualityMosaic(index).rename(index + '_post')
        # Calculate relative difference index for other indices
        rd_index = ee.Image().expression(
        '((post - pre) / sqrt((post + pre)))*100',
        {
        'pre': index_pre,
        'post': index_post,
        },).rename('rd' + index)
        
        
        #rd_index = ((index_post.subtract(index_pre)).divide(((index_post.add(index_pre))))).multiply(100).rename('rd'+ index) 
        
    return [index_pre, index_post, rd_index]


def NDVI(image):
    # Calculate Normalized Difference Vegetation Index (NDVI)
    NDVI = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
    return NDVI

def kNDVI(image):
    # Calculate Kernel NDVI (kNDVI)
    ndvi = NDVI(image)
    kNDVI = ndvi.pow(2).tanh().rename('kNDVI')
    return kNDVI

def NDWI(image):
    # Calculate Normalized Difference Water Index (NDWI)
    NDWI = image.normalizedDifference(['B3', 'B8']).rename('NDWI')
    return NDWI

def MNDWI(image):
    # Calculate Modified Normalized Difference Water Index (MNDWI)
    MNDWI = image.normalizedDifference(['B3', 'B11']).rename('MNDWI')
    return MNDWI

def RBR(image):
    # Calculate Red-Band Ratio (RBR)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    RBR = (Red.divide(Red.add(Green).add(Blue))).rename('RBR')
    return RBR

def GBR(image):
    # Calculate Green-Band Ratio (GBR)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    GBR = (Green.divide(Red.add(Green).add(Blue))).rename('GBR')
    return GBR

def EGI(image):
    # Calculate Excess Green Index (EGI)
    Red = image.select('B4')
    Green = image.select('B3')
    Blue = image.select('B2')
    EGI = (Green.multiply(ee.Number(2))).subtract(Red).subtract(Blue).rename('EGI')
    return EGI

def GRVI(image):
    # Calculate Green-Red Vegetation Index (GRVI)
    Red = image.select('B4')
    Green = image.select('B3')
    GRVI = (Green.subtract(Red)).divide(Green.add(Red)).rename('GRVI')
    return GRVI

def NDBRBI(image):
    # Calculate Sd Difference Bareness Index (NDBRBI)
    Red = image.select('B4')
    Blue = image.select('B2')
    NDBRBI = (Blue.subtract(Red)).divide(Blue.add(Red)).rename('NDBRBI')
    return NDBRBI

def SAVI(image):
    # Calculate Soil-Adjusted Vegetation Index (SAVI)
    SAVI = image.expression(
    '((NIR - RED) / (RED + NIR + L))*(1+L)',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'L': ee.Number(0.5),
    },
    ).rename('SAVI')
    return SAVI

def TSAVI(image): 
    # Calculate Transformed Soil-Adjusted Vegetation Index (TSAVI)
    TSAVI = image.expression(
    'a * ((NIR - a*RED - b) / (RED + a * NIR - a * b))',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'a': ee.Number(0.33),
        'b': ee.Number(0.1),
    },
    ).rename('TSAVI')
    return TSAVI

def MSI(image):
    # Calculate Moisture Stress Index (MSI)
    NIR = image.select('B8')
    SWIR = image.select('B11')
    MSI = SWIR.divide(NIR).rename('MSI')
    return MSI

def LSWI(image):
    # Calculate Land Surface Water Index (LSWI)
    LSWI = image.normalizedDifference(['B8', 'B11']).rename('LSWI')    
    return LSWI

def EVI(image):
    # Calculate Enhanced Vegetation Index (EVI)
    EVI = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RedEdge - 7.5 * BLUE + 1))',
    {
        'NIR': image.select('B8'),
        'RED': image.select('B4'),
        'BLUE': image.select('B2'),
        'RedEdge': image.select('B6')
    },
    ).rename('EVI')
    return EVI

def HSV(image):
    # Convert RGB image to HSV color space and extract the Hue component
    RGB = image.select(['B12','B8','B4'])
    HSV = RGB.rgbToHsv()
    H = HSV.select('hue').rename('HSV')
    return H

def HSV_1(image):
    # Convert image to HSV color space with custom calculations for H, S, and V components
    bands = ['B11', 'B8', 'B4']
    sentinel_image = image.select(bands)

    # Calculate V and S values
    V = sentinel_image.reduce(ee.Reducer.max())
    S = V.subtract(sentinel_image.reduce(ee.Reducer.min()))

    # Calculate H values using expressions
    H_R = sentinel_image.expression(
        '((nir - red) / S) * 60 + 60',  # Calculate H_R
        {'nir': sentinel_image.select('B8'), 
         'red': sentinel_image.select('B4'), 
         'S': S}
    )
    H_G = sentinel_image.expression(
        '((red - swir) / S) * 60 + 120',  # Calculate H_G
        {'red': sentinel_image.select('B4'), 
         'swir': sentinel_image.select('B11'), 
         'S': S}
    )
    H_B = sentinel_image.expression(
        '((swir - nir) / S) * 60 + 240',  # Calculate H_B
        {'swir': sentinel_image.select('B11'), 
         'nir': sentinel_image.select('B8'), 
         'S': S}
    )

    # Calculate H based on V and conditions
    H = V.where(V.eq(sentinel_image.select('B11')), H_R) \
        .where(V.eq(sentinel_image.select('B8')), H_G) \
        .where(V.eq(sentinel_image.select('B4')), H_B) \
        .where(V.eq(sentinel_image.reduce(ee.Reducer.min())), 0) \
        .rename('HSV_1')
    
    # Combine H, V, and S into an image
    HSV_1 = ee.Image([H, V, S]).rename(['HSV_1', 'V', 'S']).toUint16()
    return HSV_1
