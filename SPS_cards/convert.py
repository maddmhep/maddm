import sys 
sys.path.append('/Users/omatt/Documents/eclipse/2.3.1')

import models.check_param_card as converter

def susyhit_to_mg5(input, output, sm_input='sm_input.dat'):
    
    #start by standard conversion
    converter.convert_to_mg5card(input, output)
    # Now add the missing susyhit information
    card = converter.ParamCard(output)
    orig = converter.ParamCard(input)
    sm_base = converter.ParamCard(sm_input)

    if not card.has_param('mass', 6):
        card.add_param('mass',(6,) , value=orig['sminputs'].get(6).value)
    if not card.has_param('mass', 23):
        card.add_param('mass',(23,) , value=orig['sminputs'].get(4).value)
    if not card.has_param('mass', 15):
        card.add_param('mass',(15,) , value=orig['sminputs'].get(7).value)
    if not card.has_param('decay', 6):
        card.add_param('decay',(6,) , value=sm_base['decay'].get(6).value)
    if not card.has_param('decay', 15):
        card.add_param('decay',(15,) , value=sm_base['decay'].get(15).value)
    if not card.has_param('decay', 23):
        card.add_param('decay',(23,) , value=sm_base['decay'].get(23).value)
    if not card.has_param('decay', 24):
        card.add_param('decay',(24,) , value=sm_base['decay'].get(24).value)
    if not card.has_param('decay', 25):
        card.add_param('decay',(25,) , value=sm_base['decay'].get(25).value)
    card.write(output)

    



if __name__ == "__main__":
    input =['sps1a_slha.out','sps1b_slha.out','sps2_slha.out','sps3_slha.out','sps4_slha.out','sps5_slha.out','sps9_slha.out']

    for path in input:
        susyhit_to_mg5(path, "%s.slha2" % path)
