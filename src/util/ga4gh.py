import os
import base64
from shutil import copyfile

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_test_yml(results_payload_dict, base_xml_name, output, source, includePatientInfo, report_file):
    pmi = results_payload_dict.get("FinalReport", {}).get("PMI", {})
    sample = results_payload_dict.get("FinalReport", {}).get("Sample", {})
    variant_report = results_payload_dict.get("variant-report", {})
    biomarkers = variant_report.get("biomarkers", {})
    report_properties = results_payload_dict.get("FinalReport", {}).get("reportProperties", {})
    report_property = report_properties.get("reportProperty", [])
    genes = results_payload_dict.get("FinalReport", {}).get("Genes", [])
    properties = report_property if isinstance(report_property, list) else [report_property]
    gw_loh = next((prop for prop in properties if prop.get("@key") == "LossOfHeterozygosityScore"), None)

    os.makedirs(f"{output}/{base_xml_name}", exist_ok=True)

    yaml_file = {
        "tests": [
            {
                "name": "Foundation Medicine",
                "reference": "GRCh37",
                "sourceFile": source,
                "testType": sample.get("TestType"),
                "indexedDate": sample.get("ReceivedDate"),
                "patientIdentifier": pmi.get("MRN"),
                "patientDOB": pmi.get("DOB"),
                "patientLastName": pmi.get("LastName"),
                "bodySite": variant_report.get("@tissue-of-origin"),
                "bodySiteSystem": "http://foundation.com/bodySite",
                "bodySiteDisplay": variant_report.get("@tissue-of-origin"),
                "diagnosis": pmi.get("SubmittedDiagnosis"),
                "diagnosisDisplay": pmi.get("SubmittedDiagnosis"),
                "diagnosisSystem": "http://foundation.com/diagnosis",
                "files": [
                    {
                        "type": "shortVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{base_xml_name}/{base_xml_name}.nrm.vcf",
                    },
                    {
                        "type": "copyNumberVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{base_xml_name}/{base_xml_name}.copynumber.csv",
                    },
                    {
                        "type": "structuralVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{base_xml_name}/{base_xml_name}.structural.csv",
                    },
                ],
            }
        ]
    }

    if includePatientInfo:
        yaml_file["tests"][0]["patientInfo"] = {
            "firstName": pmi.get("FirstName"),
            "lastName": pmi.get("LastName"),
            "gender": pmi.get("Gender").lower(),
            "dob": pmi.get("DOB"),
            "identifiers": [
                {
                    "codingSystem": "http://hl7.org/fhir/v2/0203",
                    "codingCode": "MR",
                    "value": pmi.get("MRN"),
                }
            ],
        }

    if gw_loh:
        values = {
            "loh-high": "high",
            "loh-low": "low"
        }
        alterations = []
        gene = genes.get("Gene", [])
        genes_list = gene if isinstance(gene, list) else [gene]
        for gene in genes_list:
            alt = gene.get("Alterations", {}).get("Alteration", [])
            alterations.extend(alt) if isinstance(alt, list) else alterations.append(alt)

        loh_alt = next((alt for alt in alterations if values.get(alt.get("Name", "").lower())), {})

        value = gw_loh.get("value", "").replace("%", "")
        # it is possible that value is "Units Not Reported"
        yaml_file["tests"][0]["lossOfHeterozygosityScore"] = float(value) if is_number(value) else None
        yaml_file["tests"][0]["lossOfHeterozygosityStatus"] = values.get(loh_alt.get("Name", "").lower(), "unknown")

    if "microsatellite-instability" in biomarkers:
        values = {
            "MSI-H": "high",
            "MSI-L": "low",
            "MSS": "stable",
            "unknown": "indeterminate",
        }
        microsatellite_dict = biomarkers["microsatellite-instability"]
        yaml_file["tests"][0]["msi"] = values.get(
            microsatellite_dict.get("@status", "unknown")
        )

    if "tumor-mutation-burden" in biomarkers:
        tumor_dict = biomarkers["tumor-mutation-burden"]
        yaml_file["tests"][0]["tmb"] = tumor_dict.get("@status", "unknown")
        yaml_file["tests"][0]["tmbScore"] = float(tumor_dict.get("@score"))

    if report_file:
        copyfile(report_file, f"{output}/{base_xml_name}/{base_xml_name}.report.pdf")

        yaml_file["tests"][0][
            "reportFile"
        ] = f".lifeomic/foundation/{base_xml_name}/{base_xml_name}.report.pdf"

    return yaml_file
