# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 18:11:29 2021

@author: Domokos Zoltán
"""

import copy as cp

import view as vw
from equirep import Table


#Exception for command errors
class InvalidCommand(Exception):

    pass


class InvalidTableArgs(Exception):

    pass


class CmdLoop:
    
    def __init__(self, table, src_list, dst_list, dims_list, inds):
        
        self.table = table
        
        self.model_signal = self.table.base_signal(src_list, dst_list, dims_list, inds)
        
        self.view_signal = self.table.base_signal(src_list, dst_list, dims_list, inds)

        #TODO: remove this deprecated class and replace it with final sol.
        self.table_cmd_proc = DmDef() 

        return
    
    
    def set_table(self, table, src_list, dst_list, dims_list, inds):
        
        self.table = table
        
        self.model_signal = table.base_signal(src_list, dst_list, dims_list, inds)
        
        self.view_signal = table.base_signal(src_list, dst_list, dims_list, inds)
        
        return
    
    
    #TODO: remove this deprecated arg. list and replace it with final sol.
    def table_update(self, val, cmd_new, cmd_fix):
        
        self.table_cmd_proc.set_cmd(self.table, val, cmd_new, cmd_fix)
        
        return
        
    
    def set_model(self, model_signal = None, src_list = None, dst_list = None, dims_list = None, inds = None):
        
        if model_signal is None:
        
            self.model_signal = self.table.base_signal(src_list, dst_list, dims_list, inds)
            
        self.model_signal = self.table.base_change_fc(src_list, 
                                                      dst_list, 
                                                      dims_list, 
                                                      model_signal,
                                                      md = 'fc_fc')
            
        self.view_signal = self.model_signal
        
    
    
    def model_update(self, src_list, dst_list, dims_list):
        #TODO: add condition to select uniform or non-uniform versions
        self.model_signal = self.table.base_change_fc(src_list, 
                                                      dst_list, 
                                                      dims_list, 
                                                      self.model_signal,
                                                      md = 'fc_fc')
        
        self.view_signal = self.model_signal
        
        return
    
    
    def view_update(self, incr, select_bounds, did):
        
        self.view_signal.set_curr(incr, select_bounds, did)
        
        return
    
    
    def show_view(self, seldims, fix_inds, mode = 'surf'):
        
        if 'heat' == mode: 
        
            vw.plot_heatmap2(self.view_signal, seldims, fix_inds)
            
        else:
            
            vw.plot_sigfc2(self.view_signal, seldims, fix_inds)
        
        return


    def cmd_loop(self):
        
        #TODO: parse user requests and do what is needed
        
        return
    

#TODO: review and remove these deprecated classes
class DmDef:    
        
    def __init__(self):
        pass
        
    
    def cmdproc_depr(self, cmd_new, cmd_fix):
        #check the commands
        cmd = None
        cmd_set = {'t', 'f', 'n'}
        
        if (cmd_new is None) or (cmd_fix is None):
            
            raise(InvalidCommand)
        
        elif (cmd_new not in cmd_set) or (cmd_fix not in cmd_set) :
            
            raise(InvalidCommand)
        
        elif cmd_new == cmd_fix:
            
            raise(InvalidCommand)
        
        else:
                #here we have a proper command for sure
                cmd = {'src' : cmd_new, 
                       'fix' : cmd_fix}
                
                for elm in cmd_set:
                    
                    if elm not in cmd.values():
                        
                        cmd['dst'] = elm
                        
                        break 
            
        return cmd
            
    
    def set_cmd(self, table, val, cmd_new, cmd_fix):
        
        
        
        if isinstance(table, Table):
        
            newtable = cp.deepcopy(table)
            cmd = self.cmdproc_depr(cmd_new, cmd_fix)
            newtable.update(val, 0, cmd, True)
        
        else:
            #TODO: raise exception
            pass
        
        return newtable #cp.deepcopy(newtable)
    
    