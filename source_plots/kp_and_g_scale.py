#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 07:29:10 2020

@author: rodney.viereck
"""
import math

def kp_and_g_scale(HP):
	
	
#  Parameters used to Calculate Kp from Hp
	y0 = 6.9498
	A1 = 2.04984
	t1 = 1.8400

	lnval = ((HP/A1)-y0)
	if lnval <= 0.0: lnval = 1e-5
	Kp= t1*math.log(lnval)
	if Kp < 0.0: Kp = 0.

# Create G Scale from Kp

	G = 0
	if Kp>=5: G=1
	if Kp>=6: G=2
	if Kp>=7: G=3
	if Kp>=8: G=4
	if Kp>=9: G=5
	
	return(Kp, G)
