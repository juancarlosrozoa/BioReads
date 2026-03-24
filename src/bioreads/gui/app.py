"""
BioReads GUI — NGS quality control and alignment interface.
"""

from __future__ import annotations

import re
import tempfile
import threading
import webbrowser
from pathlib import Path

import customtkinter as ctk
from tkinter import filedialog
from tkinterdnd2 import TkinterDnD, DND_FILES

from bioreads.core.qc import QCEngine
from bioreads.core.aligner import AlignmentEngine
from bioreads.core.detector import detect_platform, suggest_aligner


_PLATFORMS   = ["auto", "illumina", "nanopore", "pacbio", "iontorrent", "bgi"]
_EXPERIMENTS = ["dna", "rna", "amplicon"]

_STATUS_COLORS = {
    "pending":  "#60a5fa",
    "running":  "#f59e0b",
    "done":     "#22c55e",
    "error":    "#ef4444",
}


# ── Helpers ────────────────────────────────────────────────────────────────

def _parse_drop(raw: str) -> list[Path]:
    paths = re.findall(r"\{([^}]+)\}|(\S+)", raw)
    return [Path(b or p) for b, p in paths if Path(b or p).is_file()]


# ── File row for QC list ───────────────────────────────────────────────────

class QCRow(ctk.CTkFrame):
    def __init__(self, master, path: Path, **kwargs):
        super().__init__(master, fg_color="#1e293b", corner_radius=8, **kwargs)
        self.path = path
        self.columnconfigure(1, weight=1)

        self._dot = ctk.CTkLabel(self, text="●", text_color="#60a5fa",
                                  font=ctk.CTkFont(size=14), width=24)
        self._dot.grid(row=0, column=0, padx=(10, 4), pady=6)

        self._name = ctk.CTkLabel(self, text=path.name, anchor="w",
                                   font=ctk.CTkFont(size=12, weight="bold"))
        self._name.grid(row=0, column=1, padx=4, pady=6, sticky="w")

        self._status = ctk.CTkLabel(self, text="pending",
                                     text_color="#60a5fa",
                                     font=ctk.CTkFont(size=11))
        self._status.grid(row=0, column=2, padx=(4, 12), pady=6)

        self.result = None

    def set_running(self):
        self._dot.configure(text_color="#f59e0b")
        self._status.configure(text="analyzing...", text_color="#f59e0b")

    def set_done(self, result):
        self.result = result
        self._dot.configure(text_color="#22c55e")
        q30 = f"Q30: {result.pct_q30:.1f}%"
        n50 = f"N50: {result.n50:,} bp"
        self._status.configure(
            text=f"done  {q30}  {n50}  [{result.platform}]",
            text_color="#22c55e"
        )

    def set_error(self, msg: str):
        self._dot.configure(text_color="#ef4444")
        self._status.configure(text=f"error: {msg[:60]}", text_color="#ef4444")


# ── Main application ───────────────────────────────────────────────────────

class BioReadsApp(TkinterDnD.Tk):

    def __init__(self):
        super().__init__()
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        self.title("BioReads")
        self.geometry("820x660")
        self.minsize(680, 520)

        self._qc_rows: dict[str, QCRow] = {}
        self._qc_results = []
        self._ref_path: str | None = None
        self._reads_path: str | None = None
        self._reads2_path: str | None = None

        self._build_ui()

    # ── UI ─────────────────────────────────────────────────────────────────

    def _build_ui(self):
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        # Top bar
        top = ctk.CTkFrame(self, fg_color="#0f172a", corner_radius=0)
        top.grid(row=0, column=0, sticky="ew")
        ctk.CTkLabel(top, text="BioReads",
                     font=ctk.CTkFont(size=18, weight="bold")).pack(
            side="left", padx=16, pady=12)

        # Tabs
        self._tabs = ctk.CTkTabview(self, fg_color="#0f172a",
                                     segmented_button_fg_color="#1e293b",
                                     segmented_button_selected_color="#1d4ed8",
                                     segmented_button_selected_hover_color="#1e40af")
        self._tabs.grid(row=1, column=0, sticky="nsew", padx=12, pady=(8, 0))

        self._tabs.add("  QC  ")
        self._tabs.add("  Align  ")
        self._tabs.add("  Tools  ")

        self._build_qc_tab()
        self._build_align_tab()
        self._build_tools_tab()

        # Log
        self._log = ctk.CTkTextbox(self, height=110, state="disabled",
                                    font=ctk.CTkFont(family="Courier", size=11),
                                    fg_color="#0f172a", corner_radius=0)
        self._log.grid(row=2, column=0, sticky="ew", padx=0, pady=0)

    # ── QC Tab ─────────────────────────────────────────────────────────────

    def _build_qc_tab(self):
        tab = self._tabs.tab("  QC  ")
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(1, weight=1)

        # Options row
        opts = ctk.CTkFrame(tab, fg_color="transparent")
        opts.grid(row=0, column=0, sticky="ew", pady=(0, 8))

        ctk.CTkLabel(opts, text="Platform:").pack(side="left", padx=(0, 6))
        self._qc_platform = ctk.StringVar(value="auto")
        ctk.CTkOptionMenu(opts, values=_PLATFORMS,
                          variable=self._qc_platform, width=120).pack(side="left")

        ctk.CTkLabel(opts, text="Max reads:").pack(side="left", padx=(16, 6))
        self._qc_max_reads = ctk.CTkEntry(opts, width=90)
        self._qc_max_reads.insert(0, "200000")
        self._qc_max_reads.pack(side="left")

        # Drop zone
        self._qc_drop = ctk.CTkFrame(tab, fg_color="#1e293b", corner_radius=12,
                                      border_width=2, border_color="#334155")
        self._qc_drop.grid(row=1, column=0, sticky="nsew")
        self._qc_drop.columnconfigure(0, weight=1)
        self._qc_drop.rowconfigure(0, weight=1)

        self._qc_drop_label = ctk.CTkLabel(
            self._qc_drop,
            text="Drop FASTQ files here\nor click Browse",
            text_color="#475569", font=ctk.CTkFont(size=14))
        self._qc_drop_label.grid(row=0, column=0, pady=30)

        self._qc_scroll = ctk.CTkScrollableFrame(self._qc_drop, fg_color="transparent")
        self._qc_scroll.drop_target_register(DND_FILES)
        self._qc_scroll.dnd_bind("<<Drop>>", self._qc_on_drop)

        for w in (self._qc_drop, self._qc_drop_label):
            w.drop_target_register(DND_FILES)
            w.dnd_bind("<<Drop>>", self._qc_on_drop)
            w.dnd_bind("<<DragEnter>>", lambda e: self._qc_drop.configure(border_color="#3b82f6"))
            w.dnd_bind("<<DragLeave>>", lambda e: self._qc_drop.configure(border_color="#334155"))

        # Bottom buttons
        bot = ctk.CTkFrame(tab, fg_color="#1e293b", corner_radius=0,
                            border_width=1, border_color="#334155")
        bot.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        bot.columnconfigure(2, weight=1)

        ctk.CTkButton(bot, text="Browse...", width=110, height=34,
                      command=self._qc_browse).grid(row=0, column=0, padx=(16,4), pady=10)
        ctk.CTkButton(bot, text="Clear", width=80, height=34,
                      fg_color="#334155", hover_color="#475569",
                      command=self._qc_clear).grid(row=0, column=1, padx=4, pady=10)

        self._qc_btn = ctk.CTkButton(bot, text="Run QC", width=110, height=34,
                                      command=self._run_qc)
        self._qc_btn.grid(row=0, column=3, padx=4, pady=10)

        self._qc_html_btn = ctk.CTkButton(
            bot, text="Open HTML Report", width=160, height=34,
            fg_color="#1d4ed8", hover_color="#1e40af",
            state="disabled", command=self._open_qc_html)
        self._qc_html_btn.grid(row=0, column=4, padx=(4, 16), pady=10)

        self._qc_html_path = None

    # ── Align Tab ──────────────────────────────────────────────────────────

    def _build_align_tab(self):
        tab = self._tabs.tab("  Align  ")
        tab.columnconfigure(1, weight=1)

        def row_label(text, r):
            ctk.CTkLabel(tab, text=text, anchor="e", width=120,
                         text_color="#94a3b8").grid(row=r, column=0, padx=(8,8), pady=6, sticky="e")

        # Reads
        row_label("Reads (R1 / single):", 0)
        r1_frame = ctk.CTkFrame(tab, fg_color="transparent")
        r1_frame.grid(row=0, column=1, sticky="ew", padx=(0,8), pady=6)
        r1_frame.columnconfigure(0, weight=1)
        self._r1_entry = ctk.CTkEntry(r1_frame, placeholder_text="reads.fastq")
        self._r1_entry.grid(row=0, column=0, sticky="ew")
        ctk.CTkButton(r1_frame, text="Browse", width=70, height=28,
                      command=lambda: self._pick_fastq(self._r1_entry)).grid(row=0, column=1, padx=(6,0))

        # R2
        row_label("Reads R2 (paired-end):", 1)
        r2_frame = ctk.CTkFrame(tab, fg_color="transparent")
        r2_frame.grid(row=1, column=1, sticky="ew", padx=(0,8), pady=6)
        r2_frame.columnconfigure(0, weight=1)
        self._r2_entry = ctk.CTkEntry(r2_frame, placeholder_text="reads_R2.fastq (optional)")
        self._r2_entry.grid(row=0, column=0, sticky="ew")
        ctk.CTkButton(r2_frame, text="Browse", width=70, height=28,
                      command=lambda: self._pick_fastq(self._r2_entry)).grid(row=0, column=1, padx=(6,0))

        # Reference
        row_label("Reference genome:", 2)
        ref_frame = ctk.CTkFrame(tab, fg_color="transparent")
        ref_frame.grid(row=2, column=1, sticky="ew", padx=(0,8), pady=6)
        ref_frame.columnconfigure(0, weight=1)
        self._ref_entry = ctk.CTkEntry(ref_frame, placeholder_text="genome.fa or STAR index dir")
        self._ref_entry.grid(row=0, column=0, sticky="ew")
        ctk.CTkButton(ref_frame, text="Browse", width=70, height=28,
                      command=self._pick_ref).grid(row=0, column=1, padx=(6,0))

        # Platform
        row_label("Platform:", 3)
        self._align_platform = ctk.StringVar(value="illumina")
        ctk.CTkOptionMenu(tab, values=_PLATFORMS[1:],
                          variable=self._align_platform, width=140,
                          command=self._update_aligner_label).grid(
            row=3, column=1, sticky="w", padx=(0,8), pady=6)

        # Experiment
        row_label("Experiment:", 4)
        self._align_experiment = ctk.StringVar(value="dna")
        ctk.CTkOptionMenu(tab, values=_EXPERIMENTS,
                          variable=self._align_experiment, width=140,
                          command=self._update_aligner_label).grid(
            row=4, column=1, sticky="w", padx=(0,8), pady=6)

        # Auto-selected aligner (read-only info)
        row_label("Selected aligner:", 5)
        self._aligner_label = ctk.CTkLabel(tab, text="bwa-mem2",
                                            text_color="#60a5fa",
                                            font=ctk.CTkFont(size=13))
        self._aligner_label.grid(row=5, column=1, sticky="w", padx=(0,8), pady=6)
        self._update_aligner_label()

        # Threads
        row_label("Threads:", 6)
        self._threads_entry = ctk.CTkEntry(tab, width=60)
        self._threads_entry.insert(0, "4")
        self._threads_entry.grid(row=6, column=1, sticky="w", padx=(0,8), pady=6)

        # Output
        row_label("Output file:", 7)
        out_frame = ctk.CTkFrame(tab, fg_color="transparent")
        out_frame.grid(row=7, column=1, sticky="ew", padx=(0,8), pady=6)
        out_frame.columnconfigure(0, weight=1)
        self._out_entry = ctk.CTkEntry(out_frame, placeholder_text="aligned.bam")
        self._out_entry.insert(0, "aligned.bam")
        self._out_entry.grid(row=0, column=0, sticky="ew")

        # Align button
        self._align_btn = ctk.CTkButton(tab, text="Align", height=36,
                                         command=self._run_align)
        self._align_btn.grid(row=8, column=0, columnspan=2, pady=20, padx=16, sticky="ew")

    # ── Tools Tab ──────────────────────────────────────────────────────────

    def _build_tools_tab(self):
        tab = self._tabs.tab("  Tools  ")
        tab.columnconfigure(0, weight=1)

        ctk.CTkLabel(tab, text="Installed aligners",
                     font=ctk.CTkFont(size=14, weight="bold")).pack(pady=(16, 8))

        self._tools_frame = ctk.CTkScrollableFrame(tab, fg_color="#1e293b",
                                                    corner_radius=10)
        self._tools_frame.pack(fill="both", expand=True, padx=8, pady=(0, 8))
        self._tools_frame.columnconfigure(0, weight=1)

        ctk.CTkButton(tab, text="Refresh", height=34,
                      command=self._refresh_tools).pack(pady=(0, 12))

        self._refresh_tools()

    # ── QC logic ───────────────────────────────────────────────────────────

    def _qc_on_drop(self, event):
        for p in _parse_drop(event.data):
            self._qc_add_file(p)
        self._qc_drop.configure(border_color="#334155")

    def _qc_browse(self):
        paths = filedialog.askopenfilenames(
            title="Select FASTQ files",
            filetypes=[("FASTQ", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All", "*.*")])
        for p in paths:
            self._qc_add_file(Path(p))

    def _qc_add_file(self, path: Path):
        key = str(path)
        if key in self._qc_rows:
            return
        self._qc_drop_label.grid_forget()
        self._qc_scroll.grid(row=0, column=0, sticky="nsew", padx=8, pady=8)
        row = QCRow(self._qc_scroll, path)
        row.pack(fill="x", pady=3)
        self._qc_rows[key] = row

    def _qc_clear(self):
        for r in self._qc_rows.values():
            r.destroy()
        self._qc_rows.clear()
        self._qc_results.clear()
        self._qc_html_btn.configure(state="disabled")
        self._qc_scroll.grid_forget()
        self._qc_drop_label.grid(row=0, column=0, pady=30)

    def _run_qc(self):
        if not self._qc_rows:
            self._log_write("No files added. Drop or browse FASTQ files.\n")
            return
        self._qc_btn.configure(state="disabled")
        self._qc_html_btn.configure(state="disabled")
        self._qc_results.clear()
        threading.Thread(target=self._qc_worker, daemon=True).start()

    def _qc_worker(self):
        platform = self._qc_platform.get()
        platform = None if platform == "auto" else platform
        try:
            max_reads = int(self._qc_max_reads.get())
        except ValueError:
            max_reads = 200_000

        engine = QCEngine()
        for key, row in self._qc_rows.items():
            self.after(0, row.set_running)
            try:
                result = engine.run(row.path, platform=platform, max_reads=max_reads)
                self._qc_results.append(result)
                self.after(0, row.set_done, result)
                self.after(0, self._log_write,
                           f"DONE  {row.path.name}  Q30={result.pct_q30:.1f}%  "
                           f"N50={result.n50:,}bp  [{result.platform}]\n")
            except Exception as exc:
                self.after(0, row.set_error, str(exc))
                self.after(0, self._log_write, f"ERROR {row.path.name}: {exc}\n")

        html = self._generate_qc_html()
        self._qc_html_path = html
        self.after(0, self._qc_btn.configure, {"state": "normal"})
        self.after(0, self._qc_html_btn.configure, {"state": "normal"})
        self.after(0, self._log_write, f"\nReport: {html}\n")

    def _generate_qc_html(self) -> str:
        tmp = tempfile.NamedTemporaryFile(
            suffix=".html", prefix="bioreads_qc_", delete=False)
        tmp.close()

        cards = ""
        for r in self._qc_results:
            q_color = "#22c55e" if r.pct_q30 >= 80 else "#f59e0b" if r.pct_q30 >= 60 else "#ef4444"
            cards += f"""
            <div class="card">
              <div class="card-title">{Path(r.file).name}</div>
              <div class="pill" style="background:#1e3a5f;color:#60a5fa">{r.platform}</div>
              <table class="stats">
                <tr><td>Total reads</td><td>{r.total_reads:,}</td></tr>
                <tr><td>Total bases</td><td>{r.total_bases:,}</td></tr>
                <tr><td>Length mean</td><td>{r.mean_length:.1f} bp</td></tr>
                <tr><td>N50</td><td>{r.n50:,} bp</td></tr>
                <tr><td>Mean Phred Q</td><td>{r.mean_quality:.1f}</td></tr>
                <tr><td>% bases >= Q20</td><td>{r.pct_q20:.1f}%</td></tr>
                <tr><td style="color:{q_color}">% bases >= Q30</td>
                    <td style="color:{q_color};font-weight:600">{r.pct_q30:.1f}%</td></tr>
                <tr><td>GC content</td><td>{r.gc_content:.1f}%</td></tr>
                <tr><td>Duplicate rate</td><td>{r.duplicate_rate:.1f}%</td></tr>
              </table>
              {"".join(f'<div class="adapter"><b>{k}:</b> {v} reads</div>' for k,v in r.adapter_hits.items())}
            </div>"""

        html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>BioReads QC Report</title>
<style>
  body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
       background:#0f172a;color:#e2e8f0;padding:2rem;margin:0}}
  h1{{font-size:1.5rem;margin-bottom:.25rem}}
  .sub{{color:#64748b;font-size:.9rem;margin-bottom:2rem}}
  .grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:1rem}}
  .card{{background:#1e293b;border-radius:.75rem;padding:1.25rem;border:1px solid #334155}}
  .card-title{{font-weight:700;font-size:1rem;margin-bottom:.75rem}}
  .pill{{display:inline-block;padding:.15rem .6rem;border-radius:999px;
         font-size:.75rem;margin-bottom:.75rem}}
  .stats{{width:100%;border-collapse:collapse;font-size:.85rem}}
  .stats td{{padding:.25rem .4rem}}
  .stats tr:hover td{{background:#263347}}
  .adapter{{margin-top:.5rem;font-size:.8rem;color:#f59e0b}}
</style></head>
<body>
  <h1>BioReads QC Report</h1>
  <p class="sub">{len(self._qc_results)} file(s) analyzed</p>
  <div class="grid">{cards}</div>
</body></html>"""

        Path(tmp.name).write_text(html, encoding="utf-8")
        return tmp.name

    def _open_qc_html(self):
        if self._qc_html_path:
            webbrowser.open(f"file:///{self._qc_html_path}")

    # ── Align logic ────────────────────────────────────────────────────────

    def _update_aligner_label(self, *_):
        platform = self._align_platform.get()
        experiment = self._align_experiment.get()
        info = suggest_aligner(platform, experiment)
        text = info["aligner"]
        if info.get("preset"):
            text += f"  -x {info['preset']}"
        text += f"  (level {info['level']})"
        if info.get("fallback"):
            text += f"  / fallback: {info['fallback']}"
        self._aligner_label.configure(text=text)

    def _pick_fastq(self, entry: ctk.CTkEntry):
        p = filedialog.askopenfilename(
            filetypes=[("FASTQ", "*.fastq *.fq *.fastq.gz"), ("All", "*.*")])
        if p:
            entry.delete(0, "end")
            entry.insert(0, p)

    def _pick_ref(self):
        p = filedialog.askopenfilename(
            title="Select reference FASTA",
            filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All", "*.*")])
        if p:
            self._ref_entry.delete(0, "end")
            self._ref_entry.insert(0, p)

    def _run_align(self):
        reads = self._r1_entry.get().strip()
        ref   = self._ref_entry.get().strip()
        if not reads or not ref:
            self._log_write("Reads and reference are required.\n")
            return
        self._align_btn.configure(state="disabled")
        threading.Thread(target=self._align_worker, daemon=True).start()

    def _align_worker(self):
        reads    = self._r1_entry.get().strip()
        reads2   = self._r2_entry.get().strip() or None
        ref      = self._ref_entry.get().strip()
        platform = self._align_platform.get()
        exp      = self._align_experiment.get()
        output   = self._out_entry.get().strip() or "aligned.bam"
        try:
            threads = int(self._threads_entry.get())
        except ValueError:
            threads = 4

        try:
            engine = AlignmentEngine.auto(platform=platform, experiment=exp)
            level  = engine.available_level
            if level == 0:
                info = engine.aligner_info
                self.after(0, self._log_write,
                    f"No aligner available. Install {info['aligner']} via bioconda or WSL2.\n")
                return
            self.after(0, self._log_write,
                f"Aligning with {engine.aligner_info['aligner']} (level {level})...\n")
            result = engine.align(reads, ref, output=output,
                                   reads2=reads2, threads=threads)
            self.after(0, self._log_write,
                f"Done.  Mapping rate: {result.mapping_rate:.1f}%  "
                f"({result.mapped_reads:,}/{result.total_reads:,} reads)\n"
                f"Output: {result.output_bam}\n")
        except Exception as exc:
            self.after(0, self._log_write, f"Alignment error: {exc}\n")
        finally:
            self.after(0, self._align_btn.configure, {"state": "normal"})

    # ── Tools logic ────────────────────────────────────────────────────────

    def _refresh_tools(self):
        for w in self._tools_frame.winfo_children():
            w.destroy()

        engine = AlignmentEngine("illumina")
        tools  = engine.check_installation()

        for tool, version in tools.items():
            row = ctk.CTkFrame(self._tools_frame, fg_color="#263347",
                                corner_radius=6)
            row.pack(fill="x", pady=3, padx=4)
            row.columnconfigure(1, weight=1)

            color  = "#22c55e" if version else "#475569"
            status = "installed" if version else "not found"
            ctk.CTkLabel(row, text="●", text_color=color,
                          font=ctk.CTkFont(size=13), width=24).grid(
                row=0, column=0, padx=(10, 6), pady=8)
            ctk.CTkLabel(row, text=tool, anchor="w",
                          font=ctk.CTkFont(size=12, weight="bold")).grid(
                row=0, column=1, sticky="w", pady=8)
            ctk.CTkLabel(row, text=version or status,
                          text_color=color,
                          font=ctk.CTkFont(size=11)).grid(
                row=0, column=2, padx=(4, 16), pady=8)

    # ── Log ────────────────────────────────────────────────────────────────

    def _log_write(self, text: str):
        self._log.configure(state="normal")
        self._log.insert("end", text)
        self._log.see("end")
        self._log.configure(state="disabled")


def main():
    app = BioReadsApp()
    app.mainloop()


if __name__ == "__main__":
    main()
