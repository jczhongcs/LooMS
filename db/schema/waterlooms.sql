--
-- PostgreSQL database dump
--

-- Dumped from database version 10.14 (Ubuntu 10.14-0ubuntu0.18.04.1)
-- Dumped by pg_dump version 10.14 (Ubuntu 10.14-0ubuntu0.18.04.1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: candidatepeptide; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.candidatepeptide (
    id integer NOT NULL,
    protein_id integer,
    unmodified_sequence text,
    modified_sequence text,
    charge integer,
    theoretical_mz integer,
    y_ions double precision[],
    b_ions double precision[]
);


ALTER TABLE public.candidatepeptide OWNER TO waterlooms;

--
-- Name: candidatepeptide_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.candidatepeptide_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.candidatepeptide_id_seq OWNER TO waterlooms;

--
-- Name: candidatepeptide_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.candidatepeptide_id_seq OWNED BY public.candidatepeptide.id;


--
-- Name: candidatepeptidesmodifications; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.candidatepeptidesmodifications (
    id integer NOT NULL,
    candidate_peptide_id integer,
    modification_id integer
);


ALTER TABLE public.candidatepeptidesmodifications OWNER TO waterlooms;

--
-- Name: candidatepeptidesmodifications_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.candidatepeptidesmodifications_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.candidatepeptidesmodifications_id_seq OWNER TO waterlooms;

--
-- Name: candidatepeptidesmodifications_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.candidatepeptidesmodifications_id_seq OWNED BY public.candidatepeptidesmodifications.id;


--
-- Name: candidatepeptidestrails; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.candidatepeptidestrails (
    id integer NOT NULL,
    trail_id integer,
    candidate_peptide_id integer
);


ALTER TABLE public.fragmentionstrails OWNER TO waterlooms;

--
-- Name: candidatepeptidestrails_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.candidatepeptidestrails_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.candidatepeptidestrails_id_seq OWNER TO waterlooms;

--
-- Name: candidatepeptidestrails_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.candidatepeptidestrails_id_seq OWNED BY public.candidatepeptidestrails.id;


--
-- Name: isolationwindow; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.isolationwindow (
    id integer NOT NULL,
    start_mz double precision,
    end_mz double precision,
    start_rt double precision,
    end_rt double precision
);


ALTER TABLE public.isolationwindow OWNER TO waterlooms;

--
-- Name: isolationwindow_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.isolationwindow_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.isolationwindow_id_seq OWNER TO waterlooms;

--
-- Name: isolationwindow_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.isolationwindow_id_seq OWNED BY public.isolationwindow.id;


--
-- Name: modification; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.modification (
    id integer NOT NULL,
    "position" integer,
    modification_mass double precision,
    description text
);


ALTER TABLE public.modification OWNER TO waterlooms;

--
-- Name: modification_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.modification_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.modification_id_seq OWNER TO waterlooms;

--
-- Name: modification_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.modification_id_seq OWNED BY public.modification.id;


--
-- Name: protein; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.protein (
    id integer NOT NULL,
    header text,
    protein_sequence text
);


ALTER TABLE public.protein OWNER TO waterlooms;

--
-- Name: protein_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.protein_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.protein_id_seq OWNER TO waterlooms;

--
-- Name: protein_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.protein_id_seq OWNED BY public.protein.id;


--
-- Name: trail; Type: TABLE; Schema: public; Owner: waterlooms
--

CREATE TABLE public.trail (
    id integer NOT NULL,
    isolation_window_id integer,
    local_max_mz double precision,
    local_max_rt double precision,
    intensities double precision[],
    retention_times double precision[]
);


ALTER TABLE public.trail OWNER TO waterlooms;

--
-- Name: trail_id_seq; Type: SEQUENCE; Schema: public; Owner: waterlooms
--

CREATE SEQUENCE public.trail_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.trail_id_seq OWNER TO waterlooms;

--
-- Name: trail_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: waterlooms
--

ALTER SEQUENCE public.trail_id_seq OWNED BY public.trail.id;


--
-- Name: candidatepeptide id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptide ALTER COLUMN id SET DEFAULT nextval('public.candidatepeptide_id_seq'::regclass);


--
-- Name: candidatepeptidesmodifications id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidesmodifications ALTER COLUMN id SET DEFAULT nextval('public.candidatepeptidesmodifications_id_seq'::regclass);


--
-- Name: candidatepeptidestrails id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidestrails ALTER COLUMN id SET DEFAULT nextval('public.candidatepeptidestrails_id_seq'::regclass);


--
-- Name: isolationwindow id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.isolationwindow ALTER COLUMN id SET DEFAULT nextval('public.isolationwindow_id_seq'::regclass);


--
-- Name: modification id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.modification ALTER COLUMN id SET DEFAULT nextval('public.modification_id_seq'::regclass);


--
-- Name: protein id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.protein ALTER COLUMN id SET DEFAULT nextval('public.protein_id_seq'::regclass);


--
-- Name: trail id; Type: DEFAULT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.trail ALTER COLUMN id SET DEFAULT nextval('public.trail_id_seq'::regclass);


--
-- Data for Name: candidatepeptide; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.candidatepeptide (id, protein_id, unmodified_sequence, modified_sequence, charge, theoretical_mz, y_ions, b_ions) FROM stdin;
\.


--
-- Data for Name: candidatepeptidesmodifications; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.candidatepeptidesmodifications (id, candidate_peptide_id, modification_id) FROM stdin;
\.


--
-- Data for Name: candidatepeptidestrails; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.candidatepeptidestrails (id, trail_id, candidate_peptide_id) FROM stdin;
\.


--
-- Data for Name: isolationwindow; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.isolationwindow (id, start_mz, end_mz, start_rt, end_rt) FROM stdin;
\.


--
-- Data for Name: modification; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.modification (id, "position", modification_mass, description) FROM stdin;
\.


--
-- Data for Name: protein; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.protein (id, header, protein_sequence) FROM stdin;
\.


--
-- Data for Name: trail; Type: TABLE DATA; Schema: public; Owner: waterlooms
--

COPY public.trail (id, isolation_window_id, local_max_mz, local_max_rt, intensities, retention_times) FROM stdin;
\.


--
-- Name: candidatepeptide_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.candidatepeptide_id_seq', 1, false);


--
-- Name: candidatepeptidesmodifications_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.candidatepeptidesmodifications_id_seq', 1, false);


--
-- Name: candidatepeptidestrails_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.candidatepeptidestrails_id_seq', 1, false);


--
-- Name: isolationwindow_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.isolationwindow_id_seq', 1, false);


--
-- Name: modification_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.modification_id_seq', 1, false);


--
-- Name: protein_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.protein_id_seq', 1, false);


--
-- Name: trail_id_seq; Type: SEQUENCE SET; Schema: public; Owner: waterlooms
--

SELECT pg_catalog.setval('public.trail_id_seq', 1, false);


--
-- Name: candidatepeptide candidatepeptide_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptide
    ADD CONSTRAINT candidatepeptide_pkey PRIMARY KEY (id);


--
-- Name: candidatepeptidesmodifications candidatepeptidesmodifications_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidesmodifications
    ADD CONSTRAINT candidatepeptidesmodifications_pkey PRIMARY KEY (id);


--
-- Name: candidatepeptidestrails candidatepeptidestrails_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidestrails
    ADD CONSTRAINT candidatepeptidestrails_pkey PRIMARY KEY (id);


--
-- Name: isolationwindow isolationwindow_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.isolationwindow
    ADD CONSTRAINT isolationwindow_pkey PRIMARY KEY (id);


--
-- Name: modification modification_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.modification
    ADD CONSTRAINT modification_pkey PRIMARY KEY (id);


--
-- Name: protein protein_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.protein
    ADD CONSTRAINT protein_pkey PRIMARY KEY (id);


--
-- Name: trail trail_pkey; Type: CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.trail
    ADD CONSTRAINT trail_pkey PRIMARY KEY (id);


--
-- Name: candidatepeptide_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX candidatepeptide_id_uindex ON public.candidatepeptide USING btree (id);


--
-- Name: candidatepeptidesmodifications_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX candidatepeptidesmodifications_id_uindex ON public.candidatepeptidesmodifications USING btree (id);


--
-- Name: candidatepeptidestrails_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX candidatepeptidestrails_id_uindex ON public.candidatepeptidestrails USING btree (id);


--
-- Name: isolationwindow_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX isolationwindow_id_uindex ON public.isolationwindow USING btree (id);


--
-- Name: modification_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX modification_id_uindex ON public.modification USING btree (id);


--
-- Name: protein_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX protein_id_uindex ON public.protein USING btree (id);


--
-- Name: trail_id_uindex; Type: INDEX; Schema: public; Owner: waterlooms
--

CREATE UNIQUE INDEX trail_id_uindex ON public.trail USING btree (id);


--
-- Name: candidatepeptidesmodifications candidate_peptide_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidesmodifications
    ADD CONSTRAINT candidate_peptide_id FOREIGN KEY (candidate_peptide_id) REFERENCES public.candidatepeptide(id);


--
-- Name: candidatepeptidestrails candidate_peptide_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidestrails
    ADD CONSTRAINT candidate_peptide_id FOREIGN KEY (candidate_peptide_id) REFERENCES public.candidatepeptide(id);


--
-- Name: trail isolation_window_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.trail
    ADD CONSTRAINT isolation_window_id FOREIGN KEY (isolation_window_id) REFERENCES public.isolationwindow(id);


--
-- Name: candidatepeptidesmodifications modification_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidesmodifications
    ADD CONSTRAINT modification_id FOREIGN KEY (modification_id) REFERENCES public.modification(id);


--
-- Name: candidatepeptide protein_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptide
    ADD CONSTRAINT protein_id FOREIGN KEY (protein_id) REFERENCES public.protein(id);


--
-- Name: candidatepeptidestrails trail_id; Type: FK CONSTRAINT; Schema: public; Owner: waterlooms
--

ALTER TABLE ONLY public.candidatepeptidestrails
    ADD CONSTRAINT trail_id FOREIGN KEY (trail_id) REFERENCES public.trail(id);


--
-- PostgreSQL database dump complete
--

