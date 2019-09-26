# This module contains functions that go around or fix some defects of openmc

import os
import xml.etree.ElementTree as ET

def read_periodic_surfaces():

		path_to_file = os.getcwd()+'/geometry.xml'

		tree = ET.parse(path_to_file)
		root = tree.getroot()

		periodic_surface_dict = {}

		for child in root:
			if child.tag == 'surface':
				if 'periodic_surface_id' in child.attrib:
					periodic_surface_dict[int(child.attrib['id'])] = int(child.attrib['periodic_surface_id'])

		print (periodic_surface_dict)

		return periodic_surface_dict

def add_periodic_surfaces(cell, periodic_surface_dict):

	region = cell.region
#print (region.get_surfaces())
	for surface_id in region.get_surfaces():
		if surface_id in periodic_surface_dict:
			surface = region.get_surfaces()[surface_id]
			coupled_surface_id = periodic_surface_dict[surface_id]
			coupled_surface = region.get_surfaces()[coupled_surface_id]
			surface.periodic_surface = coupled_surface
