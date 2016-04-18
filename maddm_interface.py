import madgraph.interface.master_interface as master_interface

class MadDM_interface(master_interface.MasterCmd):

    def do_helloworld(self, line):
        """print hello world"""
        print "hello world"
