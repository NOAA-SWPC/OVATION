# -*- coding: utf-8 -*-
"""
Sets the default colors for plotting


Created on Tue Mar 17 06:36:01 2020

@author: rodney.viereck
"""

def set_plot_colors():
	
	fg_color = 'yellow'	  #Set forground color
	
			
		
	land_color = (.3,.3,.2)
	ocean_color = (0., 0., .6)
	boarder_color = (.7,.7,.7)

# ****************   Set the color scale for the aurora   ***************************
#     Note, the Alpha term is for transparancy
			 
	cdict1 = {'red':  ((0.0, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.3, 0.2, 0.2),
					   (0.5, 1.0, 1.0),
					   (0.8, 1.0, 1.0),
					   (0.9, 1.0, 1.0),
				       (1.0, 0.8, 0.8)),
	
			 'green': ((0.0, 0.5, 0.5),
					   (0.1, 0.9, 0.9),
					   (0.3, 1., 1.),
					   (0.5, 1., 1.),
					   (0.8, .5, .5),
					   (0.9, 0.0, 0.0),
				       (1.0, 0.0, 0.0)),
	
			 'blue':  ((0.0, 0.0, 0.0),
						(0.1, 0.0, 0.0),
					   (0.3, 0.0, 0.0),
					   (0.5, 0.0, 0.0),
					   (0.8, 0.0, 0.0),
					   (0.9, 0.0, 0.0),
					   (1.0, 0.0, 0.0)),
						   
			'alpha': ((0.0, 0.0, .0),
	 				   (0.1, .9, .9),		
	 				   (0.3, 1.0, 1.0),
	 				   (0.5, 1.0, 1.0),
	 				   (0.8, 1.0, 1.0),
	 				   (0.9, 1.0, 1.0),
	 				   (1.0, 1.0, 1.0))
	 			}
	
	return(cdict1, fg_color, land_color, ocean_color, boarder_color)