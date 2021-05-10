import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as clr
from datetime import date
import os

def getSignificanceSymbol(pvalue, oneStar=0.05, twoStars=0.01, threeStars=0.001):
    symbol = ''
    if (pvalue<=oneStar):
        symbol = '*'
    if (pvalue<=twoStars):
        symbol = '**'
    if (pvalue<=threeStars):
        symbol = '***'
    return symbol

def createFontSizes(regular=24):
    medium = 0.8 * regular
    small = 0.7 * regular
    tiny = 0.6 * regular
    return regular, medium, small, tiny

def createTodayFormat():
    today = date.today()
    return today.strftime("%Y.%m.%d")

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def createArialNarrowFont():
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rc('font',family='Arial')
    return {'fontname':'Arial', 'stretch' : 'condensed'}

def savefigWithResolution(fig, outputDir, filePrefix, dpi=100, ext='png'):
    mkdir(outputDir)
    filename = outputDir + filePrefix + '.' + ext
    fig.savefig(filename, dpi=dpi, format=ext, bbox_inches='tight', orientation='portrait')
    print("Output " + filename)
    
def extractColorsFromColormap(n=10, colormap='jet'):
    cmap = cm.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=0, vmax=n-1) 
    return [cmap(norm(ind)) for ind in range(n)] 

def createCustomColormap(palette='white', listColors=None):
    if listColors is None:
        if palette=='black':
            listColors = ['cyan', 'royalblue', 'black', 'crimson', 'pink']
        else:
            # listColors = ['azure', 'cyan', 'black', 'crimson', 'lavenderblush']
            listColors = ['royalblue', 'cyan', 'azure', 'whitesmoke', 'lavenderblush', 'pink', 'crimson']
    return clr.LinearSegmentedColormap.from_list('custom', listColors, N=256)



