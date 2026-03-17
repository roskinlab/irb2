SEQUENCE_RECORD={
    "namespace": "roskinlab",
    "name": "sequence_record",
    "aliases": ["seq_record", "seq_rec", "record"],
    "type": "record",
    "fields": [
        {
            "name": "name",
            "type": "string"
        },
        {
            "name": "source",
            "type": ["null", "string"],
            "default": None
        },
        {
            "name": "subject",
            "type": ["null", "string"],
            "default": None
        },
        {
            "name": "sample",
            "type": ["null", "string"],
            "default": None
        },
        {
            "name": "sequence",
            "type": {
                "namespace": "roskinlab.sequence_record",
                "name": "sequence",
                "aliases": ["seq"],
                "type": "record",
                "fields": [
                    {
                        "name": "sequence",
                        "type": "string"
                    },
                    {
                        "name": "qual",
                        "type": ["string", "null"],
                        "default": None
                    },
                    {
                        "name": "annotations",
                        "type": {
                            "type": "map",
                            "values": ["null", "string", "int", "float", "boolean",
                                
                                {"namespace": "roskinlab",
                                 "name": "named_range",
                                 "type": "record",
                                 "fields": [{"name": "name",  "type": "string"},
                                            {"name": "start", "type": ["null", "int"]},
                                            {"name": "stop",  "type": ["null", "int"]}
                                           ]
                                },
                                {"namespace": "roskinlab",
                                 "name": "range",
                                 "type": "record",
                                 "fields": [{"name": "start", "type": ["null", "int"]},
                                            {"name": "stop",  "type": ["null", "int"]}
                                           ]
                                }
                            ]
                        },
                        "default": {}
                    },
                ]
            }
        },
        {
            "name": "parses",
            "default": {},
            "type": {
                "type": "map",
                "values": ["null", {
                    "namespace": "roskinlab",
                    "name": "parse",
                    "aliases": ["parse_record"],
                    "type": "record",
                    "fields": [
                        {
                            "name": "modified_sequence",
                            "type": ["null", "string"],
                            "default": None
                        },
                        {
                            "name": "chain",
                            "type": ["null", {
                                "namespace": "roskinlab.parse",
                                "name": "chain_type",
                                "aliases": ["chain"],
                                "type": "enum",
                                "symbols": ["VH", "VK", "VL", "VB", "VA", "VD", "VG"]
                            }]
                        },
                        {
                            "name": "has_stop_codon",
                            "type": ["null", "boolean"]
                        },
                        {
                            "name": "v_j_in_frame",
                            "type": ["null", "boolean"]
                        },
                        {
                            "name": "positive_strand",
                            "type": ["null", "boolean"]
                        },
                        {
                            "name": "v_frame_shift",
                            "type": ["null", "boolean"],
                            "default": None
                        },
                        {
                            "name": "alignments",
                            "type": {
                                "type": "array",
                                "items": {
                                    "namespace": "roskinlab.parse",
                                    "name": "alignment",
                                    "aliases": ["align"],
                                    "type": "record",
                                    "fields": [
                                        {
                                            "name": "type",
                                            "type": {
                                                "namespace": "roskinlab.parse",
                                                "name": "alignment_type",
                                                "aliases": ["align_type"],
                                                "type": "enum",
                                                "symbols": ["Q", "L", "V", "D", "J", "C", "T", "S"]
                                            }
                                        },
                                        {
                                            "name": "name",
                                            "type": "string"
                                        },
                                        {
                                            "name": "length",
                                            "type": "int"
                                        },
                                        {
                                            "name": "score",
                                            "type": "float"
                                        },
                                        {
                                            "name": "e_value",
                                            "type": "float"
                                        },
                                        {
                                            "name": "range",
                                            "type": "roskinlab.range"
                                        },
                                        {
                                            "name": "padding",
                                            "type": "roskinlab.range"
                                        },
                                        {
                                            "name": "alignment",
                                            "type": "string"
                                        }
                                    ]
                                }
                            }
                        },
                        {
                            "name": "annotations",
                            "type": {
                                "type": "map",
                                "values": ["null", "string", "int", "float", "boolean", "roskinlab.named_range", "roskinlab.range"]
                            },
                            "default": {}
                        }
                    ]
                }]
            }
        },
        {
            "name": "lineages",
            "type": {
                "type": "map",
                "values": "string"
            },
            "default": {}
        }
    ]
}

# pull out individual sub-records
for f in SEQUENCE_RECORD['fields']:
    if f['name'] == 'sequence':
        SEQUENCE = f['type']
    elif f['name'] == 'parses':
        PARSE = f['type']['values'][1]
assert SEQUENCE['type'] == 'record' and SEQUENCE['name'] == 'sequence'
assert PARSE['type'] == 'record' and PARSE['name'] == 'parse'
