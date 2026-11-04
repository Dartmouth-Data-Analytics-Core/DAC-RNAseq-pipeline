# Slim image for rule pca_plots (scripts/pca_plotting.py + scripts/rlog_transform.R).
# Build from the repo root:
#   docker build --platform linux/amd64 -f docker/rna_pcaplot.Dockerfile \
#     -t ghcr.io/dartmouth-data-analytics-core/rna_pcaplot:2.2 .

FROM mambaorg/micromamba:1.5.10 AS build
USER root
COPY env_config/rna_pcaplot.yaml /tmp/env.yaml
RUN micromamba create -y -q -p /opt/env -f /tmp/env.yaml \
 && micromamba clean -afy \
 # conda-forge r-base pulls a full gcc/g++/gfortran toolchain and BH headers (~1.2 GB
 # unpacked); only needed to compile R packages. Remove the ones actually installed.
 && RM="" \
 && for p in gcc_impl_linux-64 gxx_impl_linux-64 gfortran_impl_linux-64 \
             gcc_linux-64 gxx_linux-64 gfortran_linux-64 \
             binutils_impl_linux-64 binutils_linux-64 \
             sysroot_linux-64 kernel-headers_linux-64 \
             libgcc-devel_linux-64 libstdcxx-devel_linux-64 libsanitizer r-bh; do \
      ls /opt/env/conda-meta/"$p"-[0-9]*.json >/dev/null 2>&1 && RM="$RM $p"; \
    done \
 && micromamba remove -y -p /opt/env --force $RM \
 && rm -rf /opt/env/include /opt/env/man /opt/env/share/doc /opt/env/share/man \
           /opt/env/lib/R/library/*/doc /opt/env/conda-meta/*.json.bak \
 && find /opt/env -name '*.a' -delete

FROM debian:bookworm-slim
COPY --from=build /opt/env /opt/env
# PYTHONNOUSERSITE: Singularity bind-mounts $HOME, so ~/.local Python packages would shadow the image's
ENV PATH=/opt/env/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LANG=C.UTF-8 LC_ALL=C.UTF-8 MPLBACKEND=Agg PYTHONNOUSERSITE=1
CMD ["/bin/bash"]
