import os
import base64

def get_test_yml(results_payload_dict, sample_name, output, source):
    pmi = results_payload_dict.get("FinalReport", {}).get("PMI", {})
    sample = results_payload_dict.get("FinalReport", {}).get("Sample", {})
    variant_report = results_payload_dict.get("variant-report", {})
    biomarkers = variant_report.get("biomarkers", {})

    os.makedirs(f"{output}/foundation/{sample_name}", exist_ok=True)

    yaml_file = {
        "tests": [
            {
                "name": "Foundation Medicine",
                "reference": "GRCh37",
                "source": source,
                "testType": sample.get("TestType"),
                "indexedDate": sample.get("ReceivedDate"),
                "patientIdentifier": pmi.get("MRN"),
                "bodySite": variant_report.get("@tissue-of-origin"),
                "bodySiteSystem": "http://foundation.com/bodySite",
                "bodySiteDisplay": variant_report.get("@tissue-of-origin"),
                "patientInfo": {
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
                },
                "files": [
                    {
                        "type": "shortVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{sample_name}/{sample_name}.vcf",
                        "normalize": True
                    },
                    {
                        "type": "copyNumberVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{sample_name}/{sample_name}.copynumber.csv",
                    },
                    {
                        "type": "structuralVariant",
                        "sequenceType": "somatic",
                        "fileName": f".lifeomic/foundation/{sample_name}/{sample_name}.structural.csv",
                    },
                ],
            }
        ]
    }

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

    if "ReportPDF" in results_payload_dict:
        pdf = base64.b64decode(results_payload_dict["ReportPDF"])
        report_file = f"{output}/foundation/{sample_name}/{sample_name}.report.pdf"

        with open(report_file, "wb") as pdf_file:
            pdf_file.write(pdf)

        yaml_file["tests"][0][
            "reportFile"
        ] = f".lifeomic/foundation/{sample_name}/{sample_name}.report.pdf"

    return yaml_file