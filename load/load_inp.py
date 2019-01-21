#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:00:42 2018

@author: mellis
"""

import re

def rm_comment_from_line(line):
    """
    Will remove any comments in a line for parsing.
    
    Inputs:
        * line   =>  line from input file
    
    Ouputs:
        The line with comments removed
    """
    poss_comment_strs = ['#','!']
    for s in poss_comment_strs:
        words = line.split(s)
        if len(words) == 1:
            line = words[0]
        else:
            line = ''.join(words[:-1])
    return line

def str_to_num(String):
    """
    Will convert a string to a number if possible.
    
    Inputs:
        * String  =>  Any string
    
    Outpus:
        * If the string can be converted to a number it will be with ints being 
          preferable to floats. Else will return the same string
    """
    # Convert run.inp settings to bool
    if type(String) == str:
        if String.upper() == 'T':
            return True
        elif String.upper() == 'F':
            return False
    # Convert to number
    try:
        return int(String)
    except ValueError:
        try:
            return float(String)
        except ValueError:
            return String

def parse_setting_line(line):
    """
    Will parse a settings line from the run.inp file.
    
    Inputs:
        * line  =>  a line in the settings files
    
    Outputs:
        The setting, the value and any units attached to it.
    """
    words = [i for i in line.split(' ') if i]
    if len(words) == 1:
        return words[0], '', ''
    
    elif len(words) == 2:
        if '@include' in words[0]:
            return ' '.join(words), '', ''
        return words[0], str_to_num(words[1]), ''
    
    elif len(words) == 3:
        if all([j in words[1] for j in ['[',']']]):
            setting = words[0]
            unit = words[1].strip('[').strip(']')
            value = words[2]      
            return setting, str_to_num(value), unit
        else:
            return words[0], '\t'.join(words[1:]), ''
    else:
        return words[0], '\t'.join(words[1:]), ''
    

def parse_inp_file(Dict, non_nested_settings, inp_file, line_ind=0, found_sects=[]):
    """
    Will recursively parse the inp file and output this as nested dictionaries.
    
    Inputs:
        * Dict                 =>  Dict to fill with nested settings
        * non_nested_settings  =>  A dictionary of settings without any nesting
        * inp_file             =>  The inp_file to parse, this should be the 
                                    file split by '\n'. [list <str>]
        * line_ind             =>  Line to start parsing (default 0)
        * found_sects  =>  DON'T USE MAY BE DANGEROUS
    
    Ouputs:
        Each section is stored as a seperate dictionary with the nesting 
        preserved.
    """
    # Get the line and format it to remove comments etc...
    line = inp_file[line_ind]
    edit_line = line.strip()
    edit_line = rm_comment_from_line(edit_line)
    
    if edit_line or edit_line.isspace():
        # If the line is a new section or end of a section handle it
        if edit_line[0] == "&":
            edit_line = edit_line.strip("&").strip()
            if edit_line.split(' ')[0].lower() != 'end':
                section = edit_line.split(' ')[0].upper()
                found_sects.append(section)
                    
            else:
                end_sect = re.sub('end', "", edit_line, flags=re.IGNORECASE).strip().upper()
                if not end_sect:
                    end_sect = found_sects[-1]
                found_sects.remove(end_sect)
    
        
        elif found_sects:
#            print("\n")
#            print(found_sects)
#            print(edit_line)
            if Dict.get(found_sects[0]) == None:
                Dict[found_sects[0]] = {}
            curr_D = Dict[found_sects[0]]
            
            for sect in found_sects[1:]:
                if curr_D.get(sect) == None:
                    curr_D[sect] = {}
                curr_D = curr_D[sect]

            setting, value, unit = parse_setting_line(edit_line)
            curr_D[setting] = [value, unit]
            non_nested_settings[setting] = value

    if line_ind < len(inp_file)-1:
        parse_inp_file(Dict, non_nested_settings, inp_file, line_ind+1, found_sects)
    