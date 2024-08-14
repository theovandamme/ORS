import geemap  # Library for interactive mapping with Google Earth Engine
import os  # Standard library for interacting with the operating system
import ee  # Google Earth Engine library
import SAR_indices  # Custom module for SAR indices
import Geo_assets as Ga  # Custom module for geographic assets
import helper_MS as HMS  # Custom helper functions for multi-spectral processing
import MS_Indices as MSI  # Custom module for multi-spectral indices

def s2_preproc(params):
    """
    Preprocess Sentinel-2 data for the specified indices and parameters.

    Parameters:
    params (dict): Dictionary containing preprocessing parameters.

    Returns:
    list: List of images (pre-event, post-event, and relative difference index).
    """
    # Extract parameters from the dictionary
    CLOUD_CORRECTION = params['CLOUD_CORRECTION']
    START_DATE = params['START_DATE']
    EVENT_DATE = params['EVENT_DATE']
    STOP_DATE = params['STOP_DATE']
    ROI = params['ROI']
    CLIP_TO_ROI = params['CLIP_TO_ROI']
    MAKE_MAP = params['MAKE_MAP']
    SAVE_ASSET = params['SAVE_ASSET']
    INDEX = params['INDEX']
    Filename = params['FILENAME']

    ###########################################
    # 0. CHECK PARAMETERS
    ###########################################
    # Set default value for cloud correction if not specified
    if CLOUD_CORRECTION is None:
        CLOUD_CORRECTION = True

    # List of supported indices
    index_required = ['NDVI', 'kNDVI', 'NDWI', 'MNDWI', 'RBR', 'GBR', 'EGI', \
                      'GRVI', 'NDBRBI', 'SAVI', 'TSAVI', 'MSI', 'LSWI', 'EVI', 'HSV', 'HSV_1']
    # Raise an error if the specified index is not supported
    if (INDEX not in index_required):
        raise ValueError("ERROR!!! Parameter INDEX not correctly defined")

    ###########################################
    # 1. DATA SELECTION
    ###########################################
    # Select pre-event and post-event Sentinel-2 image collections
    pre = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')\
        .filterDate(START_DATE, EVENT_DATE) \
        .filterBounds(ROI)
    post = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')\
        .filterDate(EVENT_DATE, STOP_DATE) \
        .filterBounds(ROI)
    
    # Print the number of images in each collection
    print('Number of images in collection: ', pre.size().getInfo())
    print('Number of images in collection: ', post.size().getInfo()) 

    # Apply cloud correction if specified
    if CLOUD_CORRECTION:
        pre = pre.map(HMS.mask_s2_clouds)
        post = post.map(HMS.mask_s2_clouds)

    # Compute the specified index for both pre-event and post-event images
    collection = MSI.change(pre, post, INDEX)
    index_pre = collection[0]
    index_post = collection[1]
    change = collection[2]

    # Define visualization parameters for the change detection map based on the index
    palette = ['#800000', '#FF0000', '#FFA500', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF']
    if INDEX in ['NDVI', 'kNDVI', 'RBR', 'GBR', 'EGI', 'GRVI', 'NDBRBI', 'SAVI', 'TSAVI', 'MSI', 'LSWI', 'EVI']:
        visuals = {
            'bands': [change.bandNames().getInfo()[0]],
            'min': -50,
            'max': 50,
            'palette': palette
        }
    elif INDEX == 'NDWI':
        visuals = {
            'bands': [change.bandNames().getInfo()[0]],
            'min': -100,
            'max': 100,
            'palette': palette
        }
    elif INDEX == 'MNDWI':
        visuals = {
            'bands': [change.bandNames().getInfo()[0]],
            'min': -1,
            'max': 2,
            'palette': palette
        }
    elif INDEX == 'HSV_1':
        visuals = {
            'bands': [INDEX],
            'min': -100,
            'max': 100,
            'palette': palette
        }

    # Clip the images to the region of interest (ROI) if specified
    if CLIP_TO_ROI:
        clipped_collection = [i.clip(ROI) for i in collection]
        index_pre = clipped_collection[0]
        index_post = clipped_collection[1]
        change = clipped_collection[2]
    else:
        clipped_collection = collection

    # Save the processed images to Google Earth Engine assets if specified
    if SAVE_ASSET:
        names = [INDEX + '_pre', INDEX + '_post', 'rd' + INDEX]
        for number, img in enumerate(collection):
            description = Filename + '_' + names[number]
            assetId = 'projects/ee-thvdamme/assets' + '/' + Filename

            task = ee.batch.Export.image.toDrive(
                image=img,
                #assetId=assetId,
                description=description,
                region=ROI.getInfo()['coordinates'],
                scale=10,
                maxPixels=1e13
            )
            task.start()
            print(f'Exporting {description}')

    # Create a map and add the pre-event, post-event, and change detection layers
    Map = geemap.Map()
    Map.addLayer(index_pre, name=INDEX + '_pre')
    Map.addLayer(index_post, name=INDEX + '_post')
    Map.addLayer(change, vis_params=visuals, name='rd' + INDEX)

    return Map
