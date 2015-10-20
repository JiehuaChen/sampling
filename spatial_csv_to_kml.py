#!/usr/bin/python

import os, sys, glob, shutil, getopt, re
import csv
from simplekml import Kml
import pyproj


def main():				
    csv_files = glob.glob("output/*.csv")
    for f in csv_files:
        print "Writing KML file for %s" % (f,)
        csv_to_kml(f)	
        print "Done."

def csv_to_kml(input_filename):

    # open input file
    csv_file = open(input_filename,'rU')
    reader = csv.DictReader(csv_file)
    # preamble 
    input_filename_base, input_filename_ext = os.path.splitext(input_filename)
    print(input_filename_base)	
    # open output file
    kml_file = open(input_filename_base + '.kml','w')

    kml_file.write(r"""<?xml version="1.0" encoding="utf-8" ?>
                   <kml xmlns="http://www.opengis.net/kml/2.2">
                   """)
    kml_file.write("<Document><name>%s</name>" % input_filename_base)
    kml_file.write(r""" <Style id="grid1k"><IconStyle> <Icon> <color>ff0000</color> </Icon> </IconStyle></Style>""")

    kml_file.write(r"""
                   <Schema name="sample" id="sample">
                   <SimpleField name="Name" type="string"></SimpleField>
                   <SimpleField name="Description" type="string"></SimpleField>
                   <SimpleField name="GID" type="string"></SimpleField>
                   </Schema>
                   """)
    gids_unique = set()
    gids = []
    locs_1k = []        
    # main loop 
    for line in reader:

        kml_file.write('  <Placemark>\n')
        kml_file.write('  <name>GID=%s</name>\n' % (line['GID_100m']))
        kml_file.write('\t<ExtendedData><SchemaData schemaUrl=\"#sample\">\n')
        kml_file.write('  <SimpleField name="GID">%s</SimpleField>\n' % (line['GID_100m']))
        kml_file.write('\t\t</SchemaData></ExtendedData>\n')
        kml_file.write("     <Point><coordinates>%s,%s</coordinates></Point>\n" % (line['x'], line['y']))
        kml_file.write('  </Placemark>\n')

        gids_unique.add(line['GID_1k'])
        gids.append(line['GID_1k'])
        locs_1k.append([line['x_1k'], line['y_1k']])

    # epilogue
    kml_file.write('\t</Document>\n\t</kml>')
    csv_file.close()
    kml_file.close()

    gids_unique = list(gids_unique)
    locs_1k_unique = []
    for gid in gids_unique:
        locs_1k_unique.append([locs_1k[k] for k, x in enumerate(map(lambda x: x==gid, gids)) if x][0])

    path, filename = os.path.split(input_filename)

    for i, loc in enumerate(locs_1k_unique):
        print(i)
        kml=Kml()
        proj_para =  "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"
        project = pyproj.Proj(proj_para)
        loc_laea = list(project(loc[0], loc[1]))
        center_pt = kml.newpoint(name=gids_unique[i], description="1k by 1k grid", coords=[(loc[0], loc[1])])
        pol = kml.newpolygon(name="1k grid", description="A pathway in Kirstenbosch",  
                             outerboundaryis=[project(loc_laea[0]-500, loc_laea[1]+500, inverse=True), 
                                              project(loc_laea[0]+500, loc_laea[1]+500, inverse=True), 
                                              project(loc_laea[0]+500, loc_laea[1]-500, inverse=True),
                                              project(loc_laea[0]-500, loc_laea[1]-500, inverse=True), 
                                              project(loc_laea[0]-500, loc_laea[1]+500, inverse=True)],
                             innerboundaryis=[project(loc_laea[0]-500, loc_laea[1]+500, inverse=True), 
                                              project(loc_laea[0]+500, loc_laea[1]+500, inverse=True), 
                                              project(loc_laea[0]+500, loc_laea[1]-500, inverse=True),
                                              project(loc_laea[0]-500, loc_laea[1]-500, inverse=True), 
                                              project(loc_laea[0]-500, loc_laea[1]+500, inverse=True)])
        pol.polystyle.color = 'ff0000ff'
        print(gids_unique[i])
        kml.save(path+"/drone_flight_1k/"+gids_unique[i]+".kml")

if __name__ == '__main__':
    main()
