"""Code for finding, downloading, and reading static coronal hole synoptic map images from NOAA produced for WHPI project."""

import requests
from bs4 import BeautifulSoup
from sunpy.coordinates import sun 
import pandas as pd
import imageio
import os
import numpy as np

def get_chmap_files():
    """
    Get list of coronal hole map files hosted at NOAA.
    """
    url = 'https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-imagery/composites/synoptic-maps/mc-intosh/ARCHIVE_CHONLY/LEVEL3GIF'
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        links = soup.find_all('a')
        gif_files = [link['href'] for link in links if 'href' in link.attrs and link['href'].endswith('.gif')]

        # get filenames
        out_list = []
        for gif in gif_files:
            gif_splits = gif.split('_')
            crnum = gif_splits[-2]
            if len(crnum)>6:
                cr=crnum[:6]
                crext=crnum[6:]
            else:
                cr=crnum
                crext=None
            download_link = f'{url}/{gif}' 
            out_list.append([gif, cr, crext, download_link])

        # return df of filenames
        out = pd.DataFrame()
        out = pd.DataFrame(out_list, columns=['filename', 'cr_number', 'cr_ext', 'download_link'])
        return out
    else:
        print(f"Failed to retrieve the page. Status code: {response.status_code}")
        return None


def get_chmap_link(date, preferred_ext='wa'):
    """
    Given date, determine carrington rotation number, then find link for downloading associated Coronal Hole map.
    """
    crnum = sun.carrington_rotation_number(date)
    crnum = int(crnum)

    # get available files
    files_df = get_chmap_files()
    match = files_df.loc[files_df['cr_number']==f'cr{crnum}']
    if len(match)>0:
        match_ext = match.loc[match['cr_ext']==preferred_ext]
        if len(match_ext)==1:
            download_link = match_ext['download_link']
        elif len(match_ext)==0:
            download_link = match['download_link'][0]
        else:
            download_link = match_ext['download_link'][0]
        download_link = download_link.values[0]
    else:
        #no match
        download_link = None
    return download_link 
        
def load_chmap(download_link, filepath=''):
    """
    Download coronal hole map from URL and save.
    """
    filename = download_link.split('/')[-1]
    if filepath=='':
        chmap_file = f'{filename}'
    else:
        chmap_file = f'{filepath}/{filename}'

    if not os.path.exists(chmap_file):
        r = requests.get(download_link, allow_redirects=True)
        if r.status_code==200:
            open(chmap_file, 'wb').write(r.content)
    if os.path.exists(chmap_file):
        gif = imageio.mimread(chmap_file)
        return gif
    else:
        print(f'File unavailable: {chmap_file}')
        return None

def get_chmap_image(date):
    """
    Read in coronal hole map image.
    """
    download_link = get_chmap_link(date)
    if download_link is not None:
        gif = load_chmap(download_link)
    gif = np.array(gif[0])
    return gif