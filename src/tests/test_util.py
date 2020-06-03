from util.cnv import extract_copy_numbers
from util.fnv import extract_fusion_variant
from util.ga4gh import get_test_yml
import xmltodict
import os
import shutil

BASE_PATH = os.path.abspath(os.path.dirname(__file__))

def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def cleanup():
    if os.path.exists(f"{BASE_PATH}/data/foundation"):
        shutil.rmtree(f"{BASE_PATH}/data/foundation")

def test_cnv():
    cleanup()
    os.makedirs(f"{BASE_PATH}/data/foundation/SA-1612348", exist_ok=True)

    xml = read_xml(f"{BASE_PATH}/data/sample.xml")

    extract_copy_numbers(xml["rr:ResultsReport"]["rr:ResultsPayload"], "SA-1612348", f"{BASE_PATH}/data")

    csv = open(f"{BASE_PATH}/data/foundation/SA-1612348/SA-1612348.copynumber.csv")
    text = csv.read()
    csv.close()

    assert text == """sample_id,gene,copy_number,status,attributes,chromosome,start_position,end_position,interpretation
SA-1612348,CDK4,44.0,amplification,"{'number-of-exons': '7 of 7', 'status': 'amplification', 'ratio': '11.63', 'interpretation': 'known'}",chr12,58093932,58188144,Pathogenic
SA-1612348,CCND3,6.0,gain,"{'number-of-exons': '5 of 5', 'status': 'amplification', 'ratio': '2.17', 'interpretation': 'known'}",chr6,41853880,41956362,Pathogenic
SA-1612348,MYC,41.0,amplification,"{'number-of-exons': '5 of 5', 'status': 'amplification', 'ratio': '10.34', 'interpretation': 'known'}",chr8,128706589,128801451,Pathogenic
SA-1612348,PIM1,6.0,gain,"{'number-of-exons': '7 of 7', 'status': 'amplification', 'ratio': '2.14', 'interpretation': 'unknown'}",chr6,37138078,37141867,Uncertain significance
SA-1612348,RAD21,7.0,gain,"{'number-of-exons': '13 of 13', 'status': 'amplification', 'ratio': '2.69', 'interpretation': 'unknown'}",chr8,117859738,117878968,Uncertain significance
"""


def test_fnv():
    cleanup()
    os.makedirs(f"{BASE_PATH}/data/foundation/SA-1612348", exist_ok=True)

    xml = read_xml(f"{BASE_PATH}/data/sample.xml")

    extract_fusion_variant(xml["rr:ResultsReport"]["rr:ResultsPayload"], "SA-1612348", f"{BASE_PATH}/data")

    csv = open(f"{BASE_PATH}/data/foundation/SA-1612348/SA-1612348.structural.csv")
    text = csv.read()
    csv.close()

    assert text == """sample_id,gene1,gene2,effect,chromosome1,start_position1,end_position1,chromosome2,start_position2,end_position2,interpretation,sequence_type,in-frame,attributes
SA-1612348,NF1,N/A,truncation,chr17,29557687,29887856,chr6,66426718,66427149,Likely pathogenic,somatic,unknown,"{'equivocal': 'false', 'supporting-read-pairs': '83'}"
"""

def test_yml():
    cleanup()
    os.makedirs(f"{BASE_PATH}/data/foundation/SA-1612348", exist_ok=True)

    xml = read_xml(f"{BASE_PATH}/data/sample.xml")

    yml = get_test_yml(xml["rr:ResultsReport"]["rr:ResultsPayload"], "SA-1612348", f"{BASE_PATH}/data", "source")

    assert yml == {
        'tests': [
            {
                'name': 'Foundation Medicine',
                'reference': 'GRCh37',
                'source': 'source',
                'testType': 'FoundationOne Heme',
                'indexedDate': '2016-07-21',
                'patientIdentifier': '12345678',
                'bodySite': 'Bone',
                'bodySiteSystem': 'http://foundation.com/bodySite',
                'bodySiteDisplay': 'Bone',
                'patientInfo': {
                    'firstName': 'Test',
                    'lastName': 'Patient',
                    'gender': 'male',
                    'dob': '2002-12-12',
                    'identifiers': [{'codingSystem': 'http://hl7.org/fhir/v2/0203', 'codingCode': 'MR', 'value': '12345678'}]},
                'files': [
                    {
                        'type': 'shortVariant',
                        'sequenceType': 'somatic',
                        'fileName': '.lifeomic/foundation/SA-1612348/SA-1612348.vcf',
                        'normalize': True
                    },
                    {
                        'type': 'copyNumberVariant',
                        'sequenceType': 'somatic',
                        'fileName': '.lifeomic/foundation/SA-1612348/SA-1612348.copynumber.csv'
                    },
                    {
                        'type': 'structuralVariant',
                        'sequenceType': 'somatic',
                        'fileName': '.lifeomic/foundation/SA-1612348/SA-1612348.structural.csv'
                    }
                ],
                'msi': 'stable',
                'tmb': 'low',
                'tmbScore': 0.73,
                'reportFile': '.lifeomic/foundation/SA-1612348/SA-1612348.report.pdf'
            }
        ]
    }