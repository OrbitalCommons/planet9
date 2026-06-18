# Render targets for the Planet Nine manim film.
#
#   make scene FILE=scenes/preface/preface_a_foundations.py SCENE=P01Ellipse
#   make file  FILE=scenes/preface/preface_a_foundations.py     # all scenes in a file
#   make preface                                                # the whole preface
#   make data                                                   # regenerate Rust JSON
#
# Quality: Q=l (480p15, default) | m (720p30) | h (1080p60) | k (4k).

Q ?= l
ROOT := $(CURDIR)
REPO := $(abspath $(CURDIR)/..)
RUN  := PATH=$(HOME)/.local/bin:$$PATH PYTHONPATH=$(ROOT)
MANIM := $(RUN) manim render -q$(Q)

PREFACE := \
  scenes/preface/preface_a_foundations.py \
  scenes/preface/preface_b_elements.py \
  scenes/preface/preface_c_dynamics.py \
  scenes/preface/preface_d_mechanism.py \
  scenes/preface/preface_e_detection.py \
  scenes/preface/preface_f_indirect_map.py

.PHONY: scene file preface papers data clean

scene:
	$(MANIM) $(FILE) $(SCENE)

file:
	$(MANIM) -a $(FILE)

preface:
	@for f in $(PREFACE); do echo "== $$f =="; $(MANIM) -a $$f || exit 1; done

papers:
	@for f in scenes/papers/*/scene.py; do echo "== $$f =="; $(MANIM) -a $$f || exit 1; done

data:
	cd $(REPO) && cargo run -p p9-viability --offline

clean:
	rm -rf $(REPO)/media
