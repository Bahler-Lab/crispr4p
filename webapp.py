#!/usr/local/bin/python2.7
# Import modules for CGI handling
import cgi, cgitb
import os
import crispr4p.crispr4p as crp

print("Content-type: text/html\n\n")

def load_bahler_template():
    src_path = os.path.dirname(__file__)
    try:
        template_file = open(src_path + '/template/bahler_template.html')
        print template_file.read()
    except IOError as err:
        print err


# Define function to generate HTML form.
def generate_form():
    print '<H3>Please, enter gene information.</H3>\n'
    print '<TABLE BORDER = 0>\n'
    print '<FORM METHOD = post ACTION = \"webapp.py\">\n'
    print 'You can either specify by gene name,<br> \n'
    print 'Gene name: <INPUT type = text name = \"name\"><br><br><br>\n'
    print 'or specify the coordinates:<br>\n'
    print 'Chromosome: <input type = "text" name="chromosome" >\n'
    print 'Coordinates: from <INPUT type="number" name="coor_lower" min="0" max="1000000">'
    print 'to <INPUT type = number name = "coor_upper" min="0" max="1000000"><br>'
    print 'e.g. chromosome = I; Coordinates from 112000 to 115000<br><br>'
    print '<input type="submit" name="action" value="Enter"><br><br>'
    print "</FORM>\n\n"
  
   

def check_webinput(name, cr, start, stop):
    print 'You are cheking: <br>'
    if name==None:
        print 'chromosome = ', cr, '<br>'
        print 'coordinates = from', start, 'to', stop, '<br>'
    else:
        print 'name = ', name, '<br>'
    

def write_html_table(table, name):
    print '<table style="width: 100%; table-layout: fixed", boarder="0">'
    print '<tr><th align="left">', name, '</th></tr>'
    for row in table:
        print '<tr><td valign="top" align="right" bgcolor="Azure"">', row, "</td>"
        print '<td  style="width: 80%; word-wrap: break-word" bgcolor="SeaShell">', table[row], '</td></tr>'
    print '</table><br>'

  

def gRNA_report(gRNA):
    gRNA_dict = {'gRNA:':   gRNA[0],
            'Position: ':   gRNA[3],
            'Afw: '     :   gRNA[1],
            'Arv: '     :   gRNA[2]}
    write_html_table(gRNA_dict, "gRNA")


def HR_DNA_report(hr_dna):
    hr_dna_dict = {'HRfw: ': hr_dna[0],
                    'HRrv: ': hr_dna[1],
                    'Deleted DNA: ':hr_dna[2]}
    write_html_table(hr_dna_dict, "HR_DNA")

def CheckingPrimers_report(primerDesigns):
    pm = primerDesigns[0]
    pm_dict = {'Check primer left: ': pm['PRIMER_LEFT_0_SEQUENCE'],
            'Check primer right: ': pm['PRIMER_RIGHT_0_SEQUENCE'], 
            'Deleted DNA product size: ': pm['PRIMER_PAIR_0_PRODUCT_SIZE'],
            'Negative result product size: ': pm['negative_result']}
    write_html_table(pm_dict, "Primers")

def ReportPrimerDesign(name=None, cr=None, start=None, end=None):
    datapath = "data/"
    FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
    COORDINATES = datapath + 'COORDINATES.txt'
    SYNONIMS = datapath + 'SYNONIMS.txt'

    check_webinput(name, cr, start, end)
    pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS)
    try:
        ansTuple = pd.runWeb(name, cr, start, end) 
    except AssertionError as err:
        print "Error: ", err

    gRNA_report(ansTuple[0])
    HR_DNA_report(ansTuple[1])
    CheckingPrimers_report(ansTuple[2])

def webapp():
    load_bahler_template()
    form = cgi.FieldStorage()
    generate_form()
    if form.has_key("action"): 
        if form.has_key("name"):
            ReportPrimerDesign(name=str(form.getvalue("name")))
        elif form.has_key("coor_upper") and form.has_key("coor_lower") and \
                form.has_key("chromosome"):
            ReportPrimerDesign(cr=str(form.getvalue("chromosome")),
                    start=str(form.getvalue("coor_lower")),
                    end=str(form.getvalue("coor_upper")))
        else:
            print 'Error: use either gene name or coordinates.'

webapp()
