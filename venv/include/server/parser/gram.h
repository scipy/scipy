/* A Bison parser, made by GNU Bison 3.7.5.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_BASE_YY_GRAM_H_INCLUDED
# define YY_BASE_YY_GRAM_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int base_yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    IDENT = 258,                   /* IDENT  */
    UIDENT = 259,                  /* UIDENT  */
    FCONST = 260,                  /* FCONST  */
    SCONST = 261,                  /* SCONST  */
    USCONST = 262,                 /* USCONST  */
    BCONST = 263,                  /* BCONST  */
    XCONST = 264,                  /* XCONST  */
    Op = 265,                      /* Op  */
    ICONST = 266,                  /* ICONST  */
    PARAM = 267,                   /* PARAM  */
    TYPECAST = 268,                /* TYPECAST  */
    DOT_DOT = 269,                 /* DOT_DOT  */
    COLON_EQUALS = 270,            /* COLON_EQUALS  */
    EQUALS_GREATER = 271,          /* EQUALS_GREATER  */
    LESS_EQUALS = 272,             /* LESS_EQUALS  */
    GREATER_EQUALS = 273,          /* GREATER_EQUALS  */
    NOT_EQUALS = 274,              /* NOT_EQUALS  */
    ABORT_P = 275,                 /* ABORT_P  */
    ABSOLUTE_P = 276,              /* ABSOLUTE_P  */
    ACCESS = 277,                  /* ACCESS  */
    ACTION = 278,                  /* ACTION  */
    ADD_P = 279,                   /* ADD_P  */
    ADMIN = 280,                   /* ADMIN  */
    AFTER = 281,                   /* AFTER  */
    AGGREGATE = 282,               /* AGGREGATE  */
    ALL = 283,                     /* ALL  */
    ALSO = 284,                    /* ALSO  */
    ALTER = 285,                   /* ALTER  */
    ALWAYS = 286,                  /* ALWAYS  */
    ANALYSE = 287,                 /* ANALYSE  */
    ANALYZE = 288,                 /* ANALYZE  */
    AND = 289,                     /* AND  */
    ANY = 290,                     /* ANY  */
    ARRAY = 291,                   /* ARRAY  */
    AS = 292,                      /* AS  */
    ASC = 293,                     /* ASC  */
    ASENSITIVE = 294,              /* ASENSITIVE  */
    ASSERTION = 295,               /* ASSERTION  */
    ASSIGNMENT = 296,              /* ASSIGNMENT  */
    ASYMMETRIC = 297,              /* ASYMMETRIC  */
    ATOMIC = 298,                  /* ATOMIC  */
    AT = 299,                      /* AT  */
    ATTACH = 300,                  /* ATTACH  */
    ATTRIBUTE = 301,               /* ATTRIBUTE  */
    AUTHORIZATION = 302,           /* AUTHORIZATION  */
    BACKWARD = 303,                /* BACKWARD  */
    BEFORE = 304,                  /* BEFORE  */
    BEGIN_P = 305,                 /* BEGIN_P  */
    BETWEEN = 306,                 /* BETWEEN  */
    BIGINT = 307,                  /* BIGINT  */
    BINARY = 308,                  /* BINARY  */
    BIT = 309,                     /* BIT  */
    BOOLEAN_P = 310,               /* BOOLEAN_P  */
    BOTH = 311,                    /* BOTH  */
    BREADTH = 312,                 /* BREADTH  */
    BY = 313,                      /* BY  */
    CACHE = 314,                   /* CACHE  */
    CALL = 315,                    /* CALL  */
    CALLED = 316,                  /* CALLED  */
    CASCADE = 317,                 /* CASCADE  */
    CASCADED = 318,                /* CASCADED  */
    CASE = 319,                    /* CASE  */
    CAST = 320,                    /* CAST  */
    CATALOG_P = 321,               /* CATALOG_P  */
    CHAIN = 322,                   /* CHAIN  */
    CHAR_P = 323,                  /* CHAR_P  */
    CHARACTER = 324,               /* CHARACTER  */
    CHARACTERISTICS = 325,         /* CHARACTERISTICS  */
    CHECK = 326,                   /* CHECK  */
    CHECKPOINT = 327,              /* CHECKPOINT  */
    CLASS = 328,                   /* CLASS  */
    CLOSE = 329,                   /* CLOSE  */
    CLUSTER = 330,                 /* CLUSTER  */
    COALESCE = 331,                /* COALESCE  */
    COLLATE = 332,                 /* COLLATE  */
    COLLATION = 333,               /* COLLATION  */
    COLUMN = 334,                  /* COLUMN  */
    COLUMNS = 335,                 /* COLUMNS  */
    COMMENT = 336,                 /* COMMENT  */
    COMMENTS = 337,                /* COMMENTS  */
    COMMIT = 338,                  /* COMMIT  */
    COMMITTED = 339,               /* COMMITTED  */
    COMPRESSION = 340,             /* COMPRESSION  */
    CONCURRENTLY = 341,            /* CONCURRENTLY  */
    CONFIGURATION = 342,           /* CONFIGURATION  */
    CONFLICT = 343,                /* CONFLICT  */
    CONNECTION = 344,              /* CONNECTION  */
    CONSTRAINT = 345,              /* CONSTRAINT  */
    CONSTRAINTS = 346,             /* CONSTRAINTS  */
    CONTENT_P = 347,               /* CONTENT_P  */
    CONTINUE_P = 348,              /* CONTINUE_P  */
    CONVERSION_P = 349,            /* CONVERSION_P  */
    COPY = 350,                    /* COPY  */
    COST = 351,                    /* COST  */
    CREATE = 352,                  /* CREATE  */
    CROSS = 353,                   /* CROSS  */
    CSV = 354,                     /* CSV  */
    CUBE = 355,                    /* CUBE  */
    CURRENT_P = 356,               /* CURRENT_P  */
    CURRENT_CATALOG = 357,         /* CURRENT_CATALOG  */
    CURRENT_DATE = 358,            /* CURRENT_DATE  */
    CURRENT_ROLE = 359,            /* CURRENT_ROLE  */
    CURRENT_SCHEMA = 360,          /* CURRENT_SCHEMA  */
    CURRENT_TIME = 361,            /* CURRENT_TIME  */
    CURRENT_TIMESTAMP = 362,       /* CURRENT_TIMESTAMP  */
    CURRENT_USER = 363,            /* CURRENT_USER  */
    CURSOR = 364,                  /* CURSOR  */
    CYCLE = 365,                   /* CYCLE  */
    DATA_P = 366,                  /* DATA_P  */
    DATABASE = 367,                /* DATABASE  */
    DAY_P = 368,                   /* DAY_P  */
    DEALLOCATE = 369,              /* DEALLOCATE  */
    DEC = 370,                     /* DEC  */
    DECIMAL_P = 371,               /* DECIMAL_P  */
    DECLARE = 372,                 /* DECLARE  */
    DEFAULT = 373,                 /* DEFAULT  */
    DEFAULTS = 374,                /* DEFAULTS  */
    DEFERRABLE = 375,              /* DEFERRABLE  */
    DEFERRED = 376,                /* DEFERRED  */
    DEFINER = 377,                 /* DEFINER  */
    DELETE_P = 378,                /* DELETE_P  */
    DELIMITER = 379,               /* DELIMITER  */
    DELIMITERS = 380,              /* DELIMITERS  */
    DEPENDS = 381,                 /* DEPENDS  */
    DEPTH = 382,                   /* DEPTH  */
    DESC = 383,                    /* DESC  */
    DETACH = 384,                  /* DETACH  */
    DICTIONARY = 385,              /* DICTIONARY  */
    DISABLE_P = 386,               /* DISABLE_P  */
    DISCARD = 387,                 /* DISCARD  */
    DISTINCT = 388,                /* DISTINCT  */
    DO = 389,                      /* DO  */
    DOCUMENT_P = 390,              /* DOCUMENT_P  */
    DOMAIN_P = 391,                /* DOMAIN_P  */
    DOUBLE_P = 392,                /* DOUBLE_P  */
    DROP = 393,                    /* DROP  */
    EACH = 394,                    /* EACH  */
    ELSE = 395,                    /* ELSE  */
    ENABLE_P = 396,                /* ENABLE_P  */
    ENCODING = 397,                /* ENCODING  */
    ENCRYPTED = 398,               /* ENCRYPTED  */
    END_P = 399,                   /* END_P  */
    ENUM_P = 400,                  /* ENUM_P  */
    ESCAPE = 401,                  /* ESCAPE  */
    EVENT = 402,                   /* EVENT  */
    EXCEPT = 403,                  /* EXCEPT  */
    EXCLUDE = 404,                 /* EXCLUDE  */
    EXCLUDING = 405,               /* EXCLUDING  */
    EXCLUSIVE = 406,               /* EXCLUSIVE  */
    EXECUTE = 407,                 /* EXECUTE  */
    EXISTS = 408,                  /* EXISTS  */
    EXPLAIN = 409,                 /* EXPLAIN  */
    EXPRESSION = 410,              /* EXPRESSION  */
    EXTENSION = 411,               /* EXTENSION  */
    EXTERNAL = 412,                /* EXTERNAL  */
    EXTRACT = 413,                 /* EXTRACT  */
    FALSE_P = 414,                 /* FALSE_P  */
    FAMILY = 415,                  /* FAMILY  */
    FETCH = 416,                   /* FETCH  */
    FILTER = 417,                  /* FILTER  */
    FINALIZE = 418,                /* FINALIZE  */
    FIRST_P = 419,                 /* FIRST_P  */
    FLOAT_P = 420,                 /* FLOAT_P  */
    FOLLOWING = 421,               /* FOLLOWING  */
    FOR = 422,                     /* FOR  */
    FORCE = 423,                   /* FORCE  */
    FOREIGN = 424,                 /* FOREIGN  */
    FORWARD = 425,                 /* FORWARD  */
    FREEZE = 426,                  /* FREEZE  */
    FROM = 427,                    /* FROM  */
    FULL = 428,                    /* FULL  */
    FUNCTION = 429,                /* FUNCTION  */
    FUNCTIONS = 430,               /* FUNCTIONS  */
    GENERATED = 431,               /* GENERATED  */
    GLOBAL = 432,                  /* GLOBAL  */
    GRANT = 433,                   /* GRANT  */
    GRANTED = 434,                 /* GRANTED  */
    GREATEST = 435,                /* GREATEST  */
    GROUP_P = 436,                 /* GROUP_P  */
    GROUPING = 437,                /* GROUPING  */
    GROUPS = 438,                  /* GROUPS  */
    HANDLER = 439,                 /* HANDLER  */
    HAVING = 440,                  /* HAVING  */
    HEADER_P = 441,                /* HEADER_P  */
    HOLD = 442,                    /* HOLD  */
    HOUR_P = 443,                  /* HOUR_P  */
    IDENTITY_P = 444,              /* IDENTITY_P  */
    IF_P = 445,                    /* IF_P  */
    ILIKE = 446,                   /* ILIKE  */
    IMMEDIATE = 447,               /* IMMEDIATE  */
    IMMUTABLE = 448,               /* IMMUTABLE  */
    IMPLICIT_P = 449,              /* IMPLICIT_P  */
    IMPORT_P = 450,                /* IMPORT_P  */
    IN_P = 451,                    /* IN_P  */
    INCLUDE = 452,                 /* INCLUDE  */
    INCLUDING = 453,               /* INCLUDING  */
    INCREMENT = 454,               /* INCREMENT  */
    INDEX = 455,                   /* INDEX  */
    INDEXES = 456,                 /* INDEXES  */
    INHERIT = 457,                 /* INHERIT  */
    INHERITS = 458,                /* INHERITS  */
    INITIALLY = 459,               /* INITIALLY  */
    INLINE_P = 460,                /* INLINE_P  */
    INNER_P = 461,                 /* INNER_P  */
    INOUT = 462,                   /* INOUT  */
    INPUT_P = 463,                 /* INPUT_P  */
    INSENSITIVE = 464,             /* INSENSITIVE  */
    INSERT = 465,                  /* INSERT  */
    INSTEAD = 466,                 /* INSTEAD  */
    INT_P = 467,                   /* INT_P  */
    INTEGER = 468,                 /* INTEGER  */
    INTERSECT = 469,               /* INTERSECT  */
    INTERVAL = 470,                /* INTERVAL  */
    INTO = 471,                    /* INTO  */
    INVOKER = 472,                 /* INVOKER  */
    IS = 473,                      /* IS  */
    ISNULL = 474,                  /* ISNULL  */
    ISOLATION = 475,               /* ISOLATION  */
    JOIN = 476,                    /* JOIN  */
    KEY = 477,                     /* KEY  */
    LABEL = 478,                   /* LABEL  */
    LANGUAGE = 479,                /* LANGUAGE  */
    LARGE_P = 480,                 /* LARGE_P  */
    LAST_P = 481,                  /* LAST_P  */
    LATERAL_P = 482,               /* LATERAL_P  */
    LEADING = 483,                 /* LEADING  */
    LEAKPROOF = 484,               /* LEAKPROOF  */
    LEAST = 485,                   /* LEAST  */
    LEFT = 486,                    /* LEFT  */
    LEVEL = 487,                   /* LEVEL  */
    LIKE = 488,                    /* LIKE  */
    LIMIT = 489,                   /* LIMIT  */
    LISTEN = 490,                  /* LISTEN  */
    LOAD = 491,                    /* LOAD  */
    LOCAL = 492,                   /* LOCAL  */
    LOCALTIME = 493,               /* LOCALTIME  */
    LOCALTIMESTAMP = 494,          /* LOCALTIMESTAMP  */
    LOCATION = 495,                /* LOCATION  */
    LOCK_P = 496,                  /* LOCK_P  */
    LOCKED = 497,                  /* LOCKED  */
    LOGGED = 498,                  /* LOGGED  */
    MAPPING = 499,                 /* MAPPING  */
    MATCH = 500,                   /* MATCH  */
    MATCHED = 501,                 /* MATCHED  */
    MATERIALIZED = 502,            /* MATERIALIZED  */
    MAXVALUE = 503,                /* MAXVALUE  */
    MERGE = 504,                   /* MERGE  */
    METHOD = 505,                  /* METHOD  */
    MINUTE_P = 506,                /* MINUTE_P  */
    MINVALUE = 507,                /* MINVALUE  */
    MODE = 508,                    /* MODE  */
    MONTH_P = 509,                 /* MONTH_P  */
    MOVE = 510,                    /* MOVE  */
    NAME_P = 511,                  /* NAME_P  */
    NAMES = 512,                   /* NAMES  */
    NATIONAL = 513,                /* NATIONAL  */
    NATURAL = 514,                 /* NATURAL  */
    NCHAR = 515,                   /* NCHAR  */
    NEW = 516,                     /* NEW  */
    NEXT = 517,                    /* NEXT  */
    NFC = 518,                     /* NFC  */
    NFD = 519,                     /* NFD  */
    NFKC = 520,                    /* NFKC  */
    NFKD = 521,                    /* NFKD  */
    NO = 522,                      /* NO  */
    NONE = 523,                    /* NONE  */
    NORMALIZE = 524,               /* NORMALIZE  */
    NORMALIZED = 525,              /* NORMALIZED  */
    NOT = 526,                     /* NOT  */
    NOTHING = 527,                 /* NOTHING  */
    NOTIFY = 528,                  /* NOTIFY  */
    NOTNULL = 529,                 /* NOTNULL  */
    NOWAIT = 530,                  /* NOWAIT  */
    NULL_P = 531,                  /* NULL_P  */
    NULLIF = 532,                  /* NULLIF  */
    NULLS_P = 533,                 /* NULLS_P  */
    NUMERIC = 534,                 /* NUMERIC  */
    OBJECT_P = 535,                /* OBJECT_P  */
    OF = 536,                      /* OF  */
    OFF = 537,                     /* OFF  */
    OFFSET = 538,                  /* OFFSET  */
    OIDS = 539,                    /* OIDS  */
    OLD = 540,                     /* OLD  */
    ON = 541,                      /* ON  */
    ONLY = 542,                    /* ONLY  */
    OPERATOR = 543,                /* OPERATOR  */
    OPTION = 544,                  /* OPTION  */
    OPTIONS = 545,                 /* OPTIONS  */
    OR = 546,                      /* OR  */
    ORDER = 547,                   /* ORDER  */
    ORDINALITY = 548,              /* ORDINALITY  */
    OTHERS = 549,                  /* OTHERS  */
    OUT_P = 550,                   /* OUT_P  */
    OUTER_P = 551,                 /* OUTER_P  */
    OVER = 552,                    /* OVER  */
    OVERLAPS = 553,                /* OVERLAPS  */
    OVERLAY = 554,                 /* OVERLAY  */
    OVERRIDING = 555,              /* OVERRIDING  */
    OWNED = 556,                   /* OWNED  */
    OWNER = 557,                   /* OWNER  */
    PARALLEL = 558,                /* PARALLEL  */
    PARAMETER = 559,               /* PARAMETER  */
    PARSER = 560,                  /* PARSER  */
    PARTIAL = 561,                 /* PARTIAL  */
    PARTITION = 562,               /* PARTITION  */
    PASSING = 563,                 /* PASSING  */
    PASSWORD = 564,                /* PASSWORD  */
    PLACING = 565,                 /* PLACING  */
    PLANS = 566,                   /* PLANS  */
    POLICY = 567,                  /* POLICY  */
    POSITION = 568,                /* POSITION  */
    PRECEDING = 569,               /* PRECEDING  */
    PRECISION = 570,               /* PRECISION  */
    PRESERVE = 571,                /* PRESERVE  */
    PREPARE = 572,                 /* PREPARE  */
    PREPARED = 573,                /* PREPARED  */
    PRIMARY = 574,                 /* PRIMARY  */
    PRIOR = 575,                   /* PRIOR  */
    PRIVILEGES = 576,              /* PRIVILEGES  */
    PROCEDURAL = 577,              /* PROCEDURAL  */
    PROCEDURE = 578,               /* PROCEDURE  */
    PROCEDURES = 579,              /* PROCEDURES  */
    PROGRAM = 580,                 /* PROGRAM  */
    PUBLICATION = 581,             /* PUBLICATION  */
    QUOTE = 582,                   /* QUOTE  */
    RANGE = 583,                   /* RANGE  */
    READ = 584,                    /* READ  */
    REAL = 585,                    /* REAL  */
    REASSIGN = 586,                /* REASSIGN  */
    RECHECK = 587,                 /* RECHECK  */
    RECURSIVE = 588,               /* RECURSIVE  */
    REF_P = 589,                   /* REF_P  */
    REFERENCES = 590,              /* REFERENCES  */
    REFERENCING = 591,             /* REFERENCING  */
    REFRESH = 592,                 /* REFRESH  */
    REINDEX = 593,                 /* REINDEX  */
    RELATIVE_P = 594,              /* RELATIVE_P  */
    RELEASE = 595,                 /* RELEASE  */
    RENAME = 596,                  /* RENAME  */
    REPEATABLE = 597,              /* REPEATABLE  */
    REPLACE = 598,                 /* REPLACE  */
    REPLICA = 599,                 /* REPLICA  */
    RESET = 600,                   /* RESET  */
    RESTART = 601,                 /* RESTART  */
    RESTRICT = 602,                /* RESTRICT  */
    RETURN = 603,                  /* RETURN  */
    RETURNING = 604,               /* RETURNING  */
    RETURNS = 605,                 /* RETURNS  */
    REVOKE = 606,                  /* REVOKE  */
    RIGHT = 607,                   /* RIGHT  */
    ROLE = 608,                    /* ROLE  */
    ROLLBACK = 609,                /* ROLLBACK  */
    ROLLUP = 610,                  /* ROLLUP  */
    ROUTINE = 611,                 /* ROUTINE  */
    ROUTINES = 612,                /* ROUTINES  */
    ROW = 613,                     /* ROW  */
    ROWS = 614,                    /* ROWS  */
    RULE = 615,                    /* RULE  */
    SAVEPOINT = 616,               /* SAVEPOINT  */
    SCHEMA = 617,                  /* SCHEMA  */
    SCHEMAS = 618,                 /* SCHEMAS  */
    SCROLL = 619,                  /* SCROLL  */
    SEARCH = 620,                  /* SEARCH  */
    SECOND_P = 621,                /* SECOND_P  */
    SECURITY = 622,                /* SECURITY  */
    SELECT = 623,                  /* SELECT  */
    SEQUENCE = 624,                /* SEQUENCE  */
    SEQUENCES = 625,               /* SEQUENCES  */
    SERIALIZABLE = 626,            /* SERIALIZABLE  */
    SERVER = 627,                  /* SERVER  */
    SESSION = 628,                 /* SESSION  */
    SESSION_USER = 629,            /* SESSION_USER  */
    SET = 630,                     /* SET  */
    SETS = 631,                    /* SETS  */
    SETOF = 632,                   /* SETOF  */
    SHARE = 633,                   /* SHARE  */
    SHOW = 634,                    /* SHOW  */
    SIMILAR = 635,                 /* SIMILAR  */
    SIMPLE = 636,                  /* SIMPLE  */
    SKIP = 637,                    /* SKIP  */
    SMALLINT = 638,                /* SMALLINT  */
    SNAPSHOT = 639,                /* SNAPSHOT  */
    SOME = 640,                    /* SOME  */
    SQL_P = 641,                   /* SQL_P  */
    STABLE = 642,                  /* STABLE  */
    STANDALONE_P = 643,            /* STANDALONE_P  */
    START = 644,                   /* START  */
    STATEMENT = 645,               /* STATEMENT  */
    STATISTICS = 646,              /* STATISTICS  */
    STDIN = 647,                   /* STDIN  */
    STDOUT = 648,                  /* STDOUT  */
    STORAGE = 649,                 /* STORAGE  */
    STORED = 650,                  /* STORED  */
    STRICT_P = 651,                /* STRICT_P  */
    STRIP_P = 652,                 /* STRIP_P  */
    SUBSCRIPTION = 653,            /* SUBSCRIPTION  */
    SUBSTRING = 654,               /* SUBSTRING  */
    SUPPORT = 655,                 /* SUPPORT  */
    SYMMETRIC = 656,               /* SYMMETRIC  */
    SYSID = 657,                   /* SYSID  */
    SYSTEM_P = 658,                /* SYSTEM_P  */
    TABLE = 659,                   /* TABLE  */
    TABLES = 660,                  /* TABLES  */
    TABLESAMPLE = 661,             /* TABLESAMPLE  */
    TABLESPACE = 662,              /* TABLESPACE  */
    TEMP = 663,                    /* TEMP  */
    TEMPLATE = 664,                /* TEMPLATE  */
    TEMPORARY = 665,               /* TEMPORARY  */
    TEXT_P = 666,                  /* TEXT_P  */
    THEN = 667,                    /* THEN  */
    TIES = 668,                    /* TIES  */
    TIME = 669,                    /* TIME  */
    TIMESTAMP = 670,               /* TIMESTAMP  */
    TO = 671,                      /* TO  */
    TRAILING = 672,                /* TRAILING  */
    TRANSACTION = 673,             /* TRANSACTION  */
    TRANSFORM = 674,               /* TRANSFORM  */
    TREAT = 675,                   /* TREAT  */
    TRIGGER = 676,                 /* TRIGGER  */
    TRIM = 677,                    /* TRIM  */
    TRUE_P = 678,                  /* TRUE_P  */
    TRUNCATE = 679,                /* TRUNCATE  */
    TRUSTED = 680,                 /* TRUSTED  */
    TYPE_P = 681,                  /* TYPE_P  */
    TYPES_P = 682,                 /* TYPES_P  */
    UESCAPE = 683,                 /* UESCAPE  */
    UNBOUNDED = 684,               /* UNBOUNDED  */
    UNCOMMITTED = 685,             /* UNCOMMITTED  */
    UNENCRYPTED = 686,             /* UNENCRYPTED  */
    UNION = 687,                   /* UNION  */
    UNIQUE = 688,                  /* UNIQUE  */
    UNKNOWN = 689,                 /* UNKNOWN  */
    UNLISTEN = 690,                /* UNLISTEN  */
    UNLOGGED = 691,                /* UNLOGGED  */
    UNTIL = 692,                   /* UNTIL  */
    UPDATE = 693,                  /* UPDATE  */
    USER = 694,                    /* USER  */
    USING = 695,                   /* USING  */
    VACUUM = 696,                  /* VACUUM  */
    VALID = 697,                   /* VALID  */
    VALIDATE = 698,                /* VALIDATE  */
    VALIDATOR = 699,               /* VALIDATOR  */
    VALUE_P = 700,                 /* VALUE_P  */
    VALUES = 701,                  /* VALUES  */
    VARCHAR = 702,                 /* VARCHAR  */
    VARIADIC = 703,                /* VARIADIC  */
    VARYING = 704,                 /* VARYING  */
    VERBOSE = 705,                 /* VERBOSE  */
    VERSION_P = 706,               /* VERSION_P  */
    VIEW = 707,                    /* VIEW  */
    VIEWS = 708,                   /* VIEWS  */
    VOLATILE = 709,                /* VOLATILE  */
    WHEN = 710,                    /* WHEN  */
    WHERE = 711,                   /* WHERE  */
    WHITESPACE_P = 712,            /* WHITESPACE_P  */
    WINDOW = 713,                  /* WINDOW  */
    WITH = 714,                    /* WITH  */
    WITHIN = 715,                  /* WITHIN  */
    WITHOUT = 716,                 /* WITHOUT  */
    WORK = 717,                    /* WORK  */
    WRAPPER = 718,                 /* WRAPPER  */
    WRITE = 719,                   /* WRITE  */
    XML_P = 720,                   /* XML_P  */
    XMLATTRIBUTES = 721,           /* XMLATTRIBUTES  */
    XMLCONCAT = 722,               /* XMLCONCAT  */
    XMLELEMENT = 723,              /* XMLELEMENT  */
    XMLEXISTS = 724,               /* XMLEXISTS  */
    XMLFOREST = 725,               /* XMLFOREST  */
    XMLNAMESPACES = 726,           /* XMLNAMESPACES  */
    XMLPARSE = 727,                /* XMLPARSE  */
    XMLPI = 728,                   /* XMLPI  */
    XMLROOT = 729,                 /* XMLROOT  */
    XMLSERIALIZE = 730,            /* XMLSERIALIZE  */
    XMLTABLE = 731,                /* XMLTABLE  */
    YEAR_P = 732,                  /* YEAR_P  */
    YES_P = 733,                   /* YES_P  */
    ZONE = 734,                    /* ZONE  */
    NOT_LA = 735,                  /* NOT_LA  */
    NULLS_LA = 736,                /* NULLS_LA  */
    WITH_LA = 737,                 /* WITH_LA  */
    MODE_TYPE_NAME = 738,          /* MODE_TYPE_NAME  */
    MODE_PLPGSQL_EXPR = 739,       /* MODE_PLPGSQL_EXPR  */
    MODE_PLPGSQL_ASSIGN1 = 740,    /* MODE_PLPGSQL_ASSIGN1  */
    MODE_PLPGSQL_ASSIGN2 = 741,    /* MODE_PLPGSQL_ASSIGN2  */
    MODE_PLPGSQL_ASSIGN3 = 742,    /* MODE_PLPGSQL_ASSIGN3  */
    UMINUS = 743                   /* UMINUS  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 235 "gram.y"

	core_YYSTYPE core_yystype;
	/* these fields must match core_YYSTYPE: */
	int			ival;
	char	   *str;
	const char *keyword;

	char		chr;
	bool		boolean;
	JoinType	jtype;
	DropBehavior dbehavior;
	OnCommitAction oncommit;
	List	   *list;
	Node	   *node;
	ObjectType	objtype;
	TypeName   *typnam;
	FunctionParameter *fun_param;
	FunctionParameterMode fun_param_mode;
	ObjectWithArgs *objwithargs;
	DefElem	   *defelt;
	SortBy	   *sortby;
	WindowDef  *windef;
	JoinExpr   *jexpr;
	IndexElem  *ielem;
	StatsElem  *selem;
	Alias	   *alias;
	RangeVar   *range;
	IntoClause *into;
	WithClause *with;
	InferClause	*infer;
	OnConflictClause *onconflict;
	A_Indices  *aind;
	ResTarget  *target;
	struct PrivTarget *privtarget;
	AccessPriv *accesspriv;
	struct ImportQual *importqual;
	InsertStmt *istmt;
	VariableSetStmt *vsetstmt;
	PartitionElem *partelem;
	PartitionSpec *partspec;
	PartitionBoundSpec *partboundspec;
	RoleSpec   *rolespec;
	PublicationObjSpec *publicationobjectspec;
	struct SelectLimit *selectlimit;
	SetQuantifier setquantifier;
	struct GroupClause *groupclause;
	MergeWhenClause *mergewhen;
	struct KeyActions *keyactions;
	struct KeyAction *keyaction;

#line 603 "gram.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int base_yyparse (core_yyscan_t yyscanner);

#endif /* !YY_BASE_YY_GRAM_H_INCLUDED  */
