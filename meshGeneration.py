import numpy as np
from scipy.spatial import Delaunay
import math

def boundingBoxPruning(boxlimits,pt3,ptL,ptR,currColor):
    
    """
    Given a boxlimit, remove points from pt3,ptL,ptR, and currColor where the 
    3D point is outside the volume of boxlimit
    
    Paramters
    ----------
    boxlimit: list of size 6 that gives us the x,y,z bounds
        idx: value
        0: xMin
        1: xMax
        2: yMin
        3: yMax
        4: zMin
        5: zMax
    
    pt3,ptL,ptR,currColor:
        numpy arrays in where we want to remove points that are out of bounds
    """    
    delCol=[]
    #iterate through all points in pt3 and check if it is out of bounds
    for i in range(pt3.shape[1]):
        currpts=pt3[:,i]
        if(currpts[0]<boxlimits[0] or currpts[0]>boxlimits[1] or 
           currpts[1]<boxlimits[2] or currpts[1]>boxlimits[3] or 
           currpts[2]<boxlimits[4] or currpts[2]>boxlimits[5]):
            delCol.append(i)
    #delete points that are out of bounds
    pt3=np.delete(pt3, delCol, axis=1)
    ptL=np.delete(ptL, delCol, axis=1)
    ptR=np.delete(ptR, delCol, axis=1)
    currColor=np.delete(currColor, delCol, axis=1)

    
def getLength(pt1,pt2):
"""
    Get length between two 2D points using distance forumla
    
    Paramters
    ----------    
    pt1: 2D coordinate
    pt2: 2D coordinate
    
    returns
    ----------    
    dist: distance between pt1 and pt2

"""
    xDelta=(pt2[0]-pt1[0])**2
    yDelta=(pt2[1]-pt1[1])**2
    dist=math.sqrt(xDelta+yDelta)
    return dist

def checkValidTri(pt1,pt2,pt3,thresh):
"""
    Checks if all 3 sides of a triangle are shorter than thresh
    
    Paramters
    ----------    
    pt1: 2D coordinate
    pt2: 2D coordinate
    pt3: 2d coordinate
    
    thresh: float
        value that determines whether a triangle side is too long
    
    returns
    ----------    
    isValid: bool
        bool value that determines if triangle meets the thresh requirements

"""
    edge1=getLength(pt1,pt2)
    edge2=getLength(pt1,pt3)
    edge3=getLength(pt2,pt3)
    
    if(edge1>thresh or edge2>thresh or edge3>thresh):
        return False
    else:
        return True

def trianglePruning(tri,triPruningThresh):
    
"""
    removes simplices from tri where one side of the triangle > triPruningThresh
    Paramters
    ----------    
    tri: Deluaney triangle
        tri is the list of triangles that creates the mesh
    
    triPruningThresh: float
        value that determines whether a triangle side is too long
   

"""
    delRow=[]
    #iterate through all of tri's simplices
    for i in range(tri.simplices.shape[0]):
        idx=tri.simplices[i]
        pt1=tri.points[idx[0]]
        pt2=tri.points[idx[1]]
        pt3=tri.points[idx[2]]
        #add curr index of invalid triangle to delRow
        if(checkValidTri(pt1,pt2,pt3,triPruningThresh)==False):
            delRow.append(i)
    #remove all triangles where it contains an invalid triangle
    tri.simplices=np.delete(tri.simplices, delRow, axis=0)
    
def createNeighborDict(simp):
"""
    create a dict of vertex:neighbor items in order to access a vertex's neighbor in O(1) time
    Paramters
    ----------    
    simp: Deluaney tirangle simplices
        list of simplices that create a mesh

    returns
    ---------
    neighborDict: dict()
        dictionary of vertex:negibhor itesm
   

"""
    neighborDict=dict()
    numSimplex=simp.shape[0]
    #iterate through all simplices
    for i in range(numSimplex):
        simplexIndices=simp[i]
        
        #add simplex[0] to dict
        currVertex=simplexIndices[0]
        if(currVertex not in neighborDict):
            neighborDict[currVertex]=set()
        neighborDict[currVertex].add(simplexIndices[1])
        neighborDict[currVertex].add(simplexIndices[2])
        
        #add simplex[1] to dict
        currVertex=simplexIndices[1]
        if(currVertex not in neighborDict):
            neighborDict[currVertex]=set()
        neighborDict[currVertex].add(simplexIndices[0])
        neighborDict[currVertex].add(simplexIndices[2])
        
        #add simplex[2] to dict
        currVertex=simplexIndices[2]
        if(currVertex not in neighborDict):
            neighborDict[currVertex]=set()
        neighborDict[currVertex].add(simplexIndices[0])
        neighborDict[currVertex].add(simplexIndices[1])
    return neighborDict

#calculating by not including the curr index as part of the avg
def smoothMesh(nDict,pt3,numOfSmooths):
"""
    smooth mesh by taking the average coordinate values of a vertex's neighbors
    
    Paramters
    ----------    
    nDict: dict()
        dictionary that contains vertex:neighbor items
    pt3: 3XN numpy array
        array that contains the 3D points of a given vertex
    numOfSmooths: int
        number of times we smooth the mesh

"""    
    for _ in range(numOfSmooths):
        #iterate through all the vertices and its neighbors
        for vert,neighbors in nDict.items():
            count=0
            xSum=0
            ySum=0
            zSum=0
            #itearte through all neighbors to calculate the average
            for idx in neighbors:
                count+=1
                xSum+=pt3[0][idx]
                ySum+=pt3[1][idx]
                zSum+=pt3[2][idx]
            #set the curr point to be the average
            pt3[0][vert]=xSum/(count*1.0)
            pt3[1][vert]=ySum/(count*1.0)
            pt3[2][vert]=zSum/(count*1.0)
