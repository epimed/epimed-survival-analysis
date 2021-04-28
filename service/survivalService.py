import numpy as np

def getSignificanceSymbolWithHr(pvalue, hr, oneStar=0.05, twoStars=0.01, threeStars=0.001):
    symbol = ''
    if (hr<1.0):
        return symbol
    if (pvalue<=oneStar):
        symbol = '*'
    if (pvalue<=twoStars):
        symbol = '**'
    if (pvalue<=threeStars):
        symbol = '***'
    return symbol

def isSignificant(pvalue, hr):
    if (pvalue is None or hr is None):
        return False
    if (pvalue<=0.05 and hr>=1.0):
        return True
    return False


def shiftToMaxSurvivalTime(survivalTime, survivalEvent, maxSurvivalTime):
    if ((survivalTime is None) or (survivalEvent is None)):
        return np.nan, np.nan
    shiftedTime = survivalTime
    shiftedEvent = survivalEvent
    if (survivalTime>maxSurvivalTime):
        shiftedTime = maxSurvivalTime
        shiftedEvent = 0.0
    return shiftedTime, shiftedEvent   

