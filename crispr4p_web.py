#!/usr/local/bin/python2.7
# Import modules for CGI handling
import cgi, cgitb
import os
import crispr4p as crp
from flup.server.fcgi import WSGIServer

print("Content-type: text/html\n\n")

def load_bahler_template():
    src_path = os.path.dirname(__file__)
    try:
        template_file = open(src_path + '/bahler_template.html')
        print template_file.read()
    except IOError as err:
        print err


# Define function to generate HTML form.
def generate_form():
    print '<H3>Please, enter gene information.</H3>\n'
    print '<TABLE BORDER = 0>\n'
    print '<FORM METHOD = post ACTION = \"cris4p_web.py\">\n'
    print 'You can either specify by gene name,<br> \n'
    print 'Gene name: <INPUT type = text name = \"name\"><br><br><br>\n'
    print 'or specify the coordinates:<br>\n'
    print 'Chromosome: <input type = "text" name="chromosome" >\n'
    print 'Coordinates: from <INPUT type="number" name="coor_lower" min="0" max="1000000">'
    print 'to <INPUT type = number name = "coor_upper" min="0" max="1000000"><br>'
    print 'e.g. chromosome = I; Coordinates from 112000 to 115000<br><br>'
    print '<center><input type="submit" name="action" value="Enter"></center><br><br>'
    print "</FORM>\n\n"
   

def check_webinput(cr, start, stop):
    print 'You are cheking: <br>'
    print 'chromosome = ', cr, '<br>'
    print 'coordinates = from', start, 'to', stop, '<br>'
    

def gRNA_report(gRNA):
    print 'gRNA: ', gRNA[0], 'pos:', gRNA[3], "<br>"
    print 'gRNAfw: ', gRNA[1], '<br>'
    print 'gRNArv: ', gRNA[2], '<br>'


def HR_DNA_report(hr_dna):
    print 'HRfw: ', hr_dna[0], '<br>'
    print 'HRrv: ', hr_dna[1], '<br>'
    print 'Deleted DNA: ', hr_dna[2], '<br>'

def CheckingPrimers_report(primerDesigns):
    pm = primerDesigns[0]
    print 'Check primer left: <br>', pm['PRIMER_LEFT_0_SEQUENCE'], '<br>'
    print 'Check primer right: ', pm['PRIMER_RIGHT_0_SEQUENCE'], '<br>'
    print 'Deleted DNA product size: ', pm['PRIMER_PAIR_0_PRODUCT_SIZE'], '<br>'
    print 'Negative result product size: ', pm['negative_result'], '<br>'

def ReportPrimerDesign(cr, start, end):
    check_webinput(cr, start, end)
    pd = crp.PrimerDesign(crp.FASTA, crp.COORDINATES, crp.SYNONIMS)
    try:
        ansTuple = pd.runWeb(cr=cr, start=start, end=end) 
    except AssertionError as err:
        print "Error: ", err

    print '<H3>Report:</H3><br>'
    gRNA_report(ansTuple[0])
    HR_DNA_report(ansTuple[1])
    CheckingPrimers_report(ansTuple[2])

def main():
    load_bahler_template()
    form = cgi.FieldStorage()
    generate_form()
    if form.has_key("action"): 
        if form.has_key("name"):
            print 'function not ready currently\n  '
        elif form.has_key("coor_upper") and form.has_key("coor_lower") and \
                form.has_key("chromosome"):
            ReportPrimerDesign(cr=str(form.getvalue("chromosome")),
                    start=str(form.getvalue("coor_lower")),
                    end=str(form.getvalue("coor_upper")))
        else: print 'Please either use name or coordinate\n'

        
main()

#WSGIServer(app).run()

