#!/usr/bin/env python3
"""
MS-GF+ DOWN Pipeline Enhanced Modern GUI Application
High-DPI compatible, auto-detects pipeline, enhanced graphics
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import subprocess
import threading
import os
import sys
from pathlib import Path
import platform
from datetime import datetime
import queue
import time
import shutil

# -------------------- High-DPI and Graphics Configuration --------------------
class HighDPIConfig:
    @staticmethod
    def setup_high_dpi(root):
        """Configure high-DPI support for crisp rendering"""
        try:
            # Windows DPI awareness
            if sys.platform.startswith('win'):
                try:
                    from ctypes import windll
                    windll.shcore.SetProcessDpiAwareness(1)  # Per-monitor DPI aware
                except:
                    try:
                        windll.user32.SetProcessDPIAware()  # System DPI aware
                    except:
                        pass
            
            # Tkinter scaling for high-DPI displays
            root.tk.call('tk', 'scaling', 2.0)  # Increase for 600 DPI effect
            
            # For Linux/Mac, try to detect DPI
            try:
                dpi = root.winfo_fpixels('1i')
                if dpi > 120:  # High DPI detected
                    scale_factor = min(dpi / 96.0, 2.5)  # Cap at 2.5x
                    root.tk.call('tk', 'scaling', scale_factor)
            except:
                pass
                
        except Exception as e:
            print(f"DPI setup warning: {e}")

# -------------------- Enhanced Visual Style Configuration --------------------
class ModernStyle:
    COLORS = {
        "primary": "#2563eb",        # Professional blue
        "primary_hover": "#1d4ed8",  # Darker blue
        "secondary": "#06b6d4",      # Cyan accent
        "success": "#10b981",        # Green
        "error": "#ef4444",          # Red
        "error_hover": "#dc2626",    # Darker red
        "warning": "#f59e0b",        # Amber
        "text_primary": "#0f172a",   # Very dark slate
        "text_secondary": "#475569", # Medium slate
        "bg_light": "#f8fafc",       # Very light slate
        "bg_card": "#ffffff",        # Pure white
        "bg_dark": "#1e293b",        # Dark slate for terminal
        "border": "#e2e8f0",         # Light border
        "accent": "#8b5cf6",         # Purple accent
    }

    FONTS = {
        'title': ('SF Pro Display', 24, 'bold'),     # macOS style
        'subtitle': ('SF Pro Display', 16, 'bold'), 
        'body': ('SF Pro Text', 11),                # Better readability
        'button': ('SF Pro Text', 11, 'bold'),
        'code': ('SF Mono', 10),                    # Monospace
        'small': ('SF Pro Text', 9),
    }
    
    # Fallback fonts for cross-platform
    FONT_FALLBACKS = {
        'SF Pro Display': ['SF Pro Display', 'Segoe UI', 'Helvetica Neue', 'Arial'],
        'SF Pro Text': ['SF Pro Text', 'Segoe UI', 'Helvetica Neue', 'Arial'],
        'SF Mono': ['SF Mono', 'Consolas', 'Monaco', 'Courier New'],
    }
    
    @staticmethod
    def get_font(font_key):
        """Get font with fallback support"""
        if font_key not in ModernStyle.FONTS:
            return ('Arial', 10)
            
        font_spec = ModernStyle.FONTS[font_key]
        font_name = font_spec[0]
        
        if font_name in ModernStyle.FONT_FALLBACKS:
            for fallback in ModernStyle.FONT_FALLBACKS[font_name]:
                try:
                    # Test if font exists by creating a temp font object
                    test_font = (fallback, font_spec[1], font_spec[2] if len(font_spec) > 2 else 'normal')
                    return test_font
                except:
                    continue
        
        return font_spec

# -------------------- Pipeline Detection --------------------
class PipelineDetector:
    @staticmethod
    def find_pipeline_script():
        """Auto-detect msgpypulse pipeline script location"""
        script_names = ['msgpypulse_cli.py', 'msgfpdown_pipeline.py', 'msgpypulse.py']
        
        search_paths = [
            # Current directory
            Path.cwd(),
            # Script directory
            Path(__file__).parent,
            # Common installation paths
            Path.home() / 'bin',
            Path('/usr/local/bin'),
            Path('/opt/msgfpdown'),
        ]
        
        # Add conda environment paths if available
        if 'CONDA_PREFIX' in os.environ:
            search_paths.extend([
                Path(os.environ['CONDA_PREFIX']) / 'bin',
                Path(os.environ['CONDA_PREFIX']) / 'Scripts',  # Windows conda
            ])
        
        # Check if it's available as a Python module
        # for script_name in script_names:
        #     module_name = script_name.replace('.py', '')
        #     try:
        #         result = subprocess.run([sys.executable, '-c', f'import {module_name}'], 
        #                               capture_output=True)
        #         if result.returncode == 0:
        #             return f"python -m {module_name}"
            # except:
            #     pass
        
        # Search in paths
        for path in search_paths:
            for script_name in script_names:
                full_path = path / script_name
                if full_path.is_file():
                    return str(full_path)
        
        # Check if it's in PATH
        for script_name in script_names:
            if shutil.which(script_name):
                return shutil.which(script_name)
        
        return None

# -------------------- Enhanced Tooltip --------------------
class EnhancedToolTip:
    def __init__(self, widget, text, delay=500):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.tip = None
        self.timer = None
        widget.bind("<Enter>", self.on_enter)
        widget.bind("<Leave>", self.on_leave)

    def on_enter(self, event=None):
        self.cancel_timer()
        self.timer = self.widget.after(self.delay, self.show_tip)

    def on_leave(self, event=None):
        self.cancel_timer()
        self.hide_tip()

    def cancel_timer(self):
        if self.timer:
            self.widget.after_cancel(self.timer)
            self.timer = None

    def show_tip(self):
        if self.tip:
            return
        
        x = self.widget.winfo_rootx() + 25
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 5
        
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{x}+{y}")
        self.tip.attributes('-topmost', True)
        
        # Enhanced styling
        frame = tk.Frame(self.tip, background="#2d3748", relief="solid", borderwidth=1)
        frame.pack()
        
        label = tk.Label(frame, text=self.text, justify=tk.LEFT,
                        background="#2d3748", foreground="#ffffff",
                        font=ModernStyle.get_font('small'),
                        padx=8, pady=6, borderwidth=0)
        label.pack()

    def hide_tip(self):
        if self.tip:
            self.tip.destroy()
            self.tip = None

# -------------------- Main Enhanced GUI Class --------------------
class MSGFPDOWNEnhancedGUI:
    def __init__(self, root):
        self.root = root
        self.pipeline_path = None
        self.scaling_factor = self.detect_scaling()
        self.setup_window()
        self.setup_styles()
        self.detect_pipeline()

        # Variables for required parameters
        self.sic_input_dir = tk.StringVar()
        self.final_output = tk.StringVar()
        self.database_path = tk.StringVar()

        # Variables for optional parameters (updated to match pipeline)
        self.score_field = tk.StringVar(value="MSMSScore")
        self.threshold = tk.StringVar(value="10")
        self.score_field2 = tk.StringVar()
        self.threshold2 = tk.StringVar()
        self.num_pep = tk.StringVar(value="2")
        self.mode = tk.StringVar(value="all_matches")
        self.rollup = tk.StringVar(value="sum")
        self.outlier_alpha = tk.StringVar()
        self.coverage_tsv = tk.StringVar()
        self.cluster_crosstab_tsv = tk.StringVar()
        self.run_cluster_step = tk.BooleanVar(value=False)
        
        # New parameter from pipeline
        #self.filter_method = tk.StringVar()

        self.is_running = False
        self.process = None
        self._stop_requested = False
        self._log_queue = queue.Queue()

        # Build UI
        self.create_widgets()
        self.bind_shortcuts()
        # Start periodic log updater
        self.root.after(200, self._poll_log_queue)
    
    def detect_scaling(self):
        """Detect platform and set appropriate scaling factor"""
        scaling = 1.0
        system = platform.system()
        
        try:
            if system == "Windows":
                # Enable high-DPI awareness
                from ctypes import windll
                windll.shcore.SetProcessDpiAwareness(1)  # Windows 8.1+
                scaling = 1.25  # adjust based on monitor
            elif system == "Darwin":
                # macOS usually handles scaling automatically
                scaling = 1.0
            else:
                # Linux / others
                scaling = 1.25  # tweak as needed
        except Exception:
            pass
        
        return scaling

    def detect_pipeline(self):
        """Detect and validate pipeline script"""
        self.pipeline_path = PipelineDetector.find_pipeline_script()
        if not self.pipeline_path:
            messagebox.showerror(
                "Pipeline Not Found", 
                "Could not locate mspypulse.py or msgpypulse.py script.\n\n"
                "Please ensure the pipeline script is:\n"
                "• In the current directory\n"
                "• In your conda environment bin\n"
                "• Installed as a Python module\n"
                "• Available in your PATH"
            )
            self.root.quit()

    def setup_window(self):
        """Configure main window with DPI scaling"""
        self.root.title("MSGFPDOWN Pipeline - Enhanced GUI")
        width, height = int(1400 * self.scaling_factor), int(1000 * self.scaling_factor)
        min_width, min_height = int(1200 * self.scaling_factor), int(800 * self.scaling_factor)
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(min_width, min_height)
        self.root.configure(bg="#f0f0f0")  # Example bg color

        # Apply Tk scaling
        self.root.tk.call('tk', 'scaling', self.scaling_factor)

        # Center window
        self.root.update_idletasks()
        x = (self.root.winfo_screenwidth() // 2) - (self.root.winfo_width() // 2)
        y = (self.root.winfo_screenheight() // 2) - (self.root.winfo_height() // 2)
        self.root.geometry(f"+{x}+{y}")

    def setup_styles(self):
        """Enhanced styling with DPI-aware fonts and ttk themes"""
        self.style = ttk.Style()
        
        # Pick best available theme
        available_themes = self.style.theme_names()
        preferred_themes = ["clam", "alt", "aqua"]
        for theme in preferred_themes:
            if theme in available_themes:
                self.style.theme_use(theme)
                break

        # Base font sizes (scaled)
        base_title = int(20 * self.scaling_factor)
        base_subtitle = int(16 * self.scaling_factor)
        base_body = int(12 * self.scaling_factor)
        base_button = int(12 * self.scaling_factor)
        font_family = "Arial"

         # Fonts
        self.title_font = (font_family, base_title, "bold")
        self.subtitle_font = (font_family, base_subtitle)
        self.body_font = (font_family, base_body)
        self.button_font = (font_family, base_button, "bold")

        # Label styles
        self.style.configure("Title.TLabel", font=self.title_font, foreground="#2a2a2a", background="#f0f0f0")
        self.style.configure("Subtitle.TLabel", font=self.subtitle_font, foreground="#333333", background="#f0f0f0")
        self.style.configure("Body.TLabel", font=self.body_font, foreground="#222222", background="#ffffff")

        # Frame styles
        self.style.configure("Card.TFrame",
            background="#ffffff",
            relief="flat",
            borderwidth=1
        )

        self.style.configure("Header.TFrame",
            background="#007acc",
            relief="flat"
        )

        # Button styles
        self.style.configure("Primary.TButton",
            font=self.button_font,
            foreground="white",
            background="#007acc",
            borderwidth=0
        )
        self.style.map("Primary.TButton",
            background=[("!disabled", "#007acc"), ("active", "#005f99"), ("disabled", "#cccccc")],
            foreground=[("!disabled", "white"), ("disabled", "#888888")]
        )

        self.style.configure("Success.TButton",
            font=self.button_font,
            foreground="white",
            background="#28a745",  # green
            borderwidth=0
        )
        self.style.map("Success.TButton",
            background=[("!disabled", "#28a745"), ("active", "#218838"), ("disabled", "#cccccc")],
            foreground=[("!disabled", "white"), ("disabled", "#888888")]
        )

        self.style.configure("Danger.TButton",
            font=("Arial", base_button, "bold"),
            foreground="white",
            background="#cc0000",
            borderwidth=0
        )
        self.style.map("Danger.TButton",
            background=[("!disabled", "#cc0000"), ("active", "#990000"), ("disabled", "#cccccc")],
            foreground=[("!disabled", "white"), ("disabled", "#888888")]
        )

        # Input fields
        self.style.configure("Modern.TEntry", font=self.body_font, fieldbackground="#ffffff", borderwidth=1)
        self.style.configure("Modern.TCombobox", font=self.body_font, fieldbackground="#ffffff", borderwidth=1)

        # Notebook
        self.style.configure("TNotebook", background="#f0f0f0", borderwidth=0)
        self.style.configure("TNotebook.Tab", font=self.button_font, padding=[20, 10])

        # New Styles 
        self.style.configure("Large.TLabelframe.Label", font=self.subtitle_font)
        self.style.configure("Large.TLabelframe", labelanchor='n')
        self.style.configure("Header.TFrame", background=ModernStyle.COLORS['primary'])
        self.style.configure("HeaderSubtitle.TLabel",
                     background=ModernStyle.COLORS['primary'],
                     foreground="white",
                     font=ModernStyle.get_font('body'))

    def create_widgets(self):
        """Create enhanced UI with better organization"""
        # Main container with padding
        main_container = ttk.Frame(self.root, padding=20)
        main_container.pack(fill=tk.BOTH, expand=True)

        # Enhanced header
        self.create_header(main_container)

        # Pipeline status indicator
        self.create_status_indicator(main_container)

        # Main content area
        self.notebook = ttk.Notebook(main_container)
        self.notebook.pack(fill=tk.BOTH, expand=True, pady=(20, 0))

        # Enhanced tabs
        self.create_configuration_tab()
        self.create_advanced_tab()
        self.create_monitoring_tab()
        self.create_help_tab()

        # Enhanced footer with better button layout
        self.create_footer(main_container)

        # Bind close event
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

    def create_header(self, parent):
        """Create enhanced header with gradient-like effect"""
        header_container = ttk.Frame(parent, height=100, style="Header.TFrame")
        header_container.pack(fill=tk.X, pady=(0, 20))
        header_container.pack_propagate(False)

        # Main header frame with gradient simulation
        header_frame = ttk.Frame(header_container,  style="Header.TFrame")
        header_frame.pack(fill=tk.BOTH, expand=True)

        # Content frame
        content_frame = ttk.Frame(header_frame,  style="Header.TFrame")
        content_frame.pack(fill=tk.BOTH, expand=True, padx=30, pady=20)

        # Left side - Logo and text
        left_frame = ttk.Frame(content_frame,  style="Header.TFrame")
        left_frame.pack(side=tk.LEFT, fill=tk.Y)

        # Enhanced logo
        logo_frame = ttk.Frame(left_frame,  style="Header.TFrame")
        logo_frame.pack(side=tk.LEFT, padx=(0, 20))
        
        protein_canvas = tk.Canvas(
            logo_frame,
            width=60,
            height=60,
            bg=ModernStyle.COLORS['primary'],
            highlightthickness=0
        )
        protein_canvas.pack()

        nodes = [(20, 25), (40, 15), (60, 25), (40, 40), (25, 55), (55, 55)]
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                x1, y1 = nodes[i]
                x2, y2 = nodes[j]
                protein_canvas.create_line(x1, y1, x2, y2, fill="white", width=2)

        for x, y in nodes:
            protein_canvas.create_oval(x-4, y-4, x+4, y+4, fill="white", outline="white")

        # Title and subtitle
        text_frame = ttk.Frame(left_frame, style="Header.TFrame")
        text_frame.pack(side=tk.LEFT, fill=tk.Y, expand=True, padx=(10,0))
        
        title_label = ttk.Label(text_frame, text="MS-GF+ Downstream Pipeline", 
                              font=ModernStyle.get_font('title'),
                              foreground='white', background=ModernStyle.COLORS['primary'])
        title_label.pack(anchor='w')
        
        subtitle_label = ttk.Label(text_frame, text="Proteomics Data Analysis Interface",
                           style="HeaderSubtitle.TLabel")
        subtitle_label.pack(anchor='w')

        # Right side - Pipeline info
        right_frame = tk.Frame(content_frame, bg=ModernStyle.COLORS['primary'])
        right_frame.pack(side=tk.RIGHT, fill=tk.Y)
        
        pipeline_info = tk.Label(right_frame, 
                               text=f"Pipeline: {os.path.basename(self.pipeline_path) if self.pipeline_path else 'Not Found'}", 
                               font=ModernStyle.get_font('small'),
                               fg='white', bg=ModernStyle.COLORS['primary'])
        pipeline_info.pack(anchor='e')

    def create_status_indicator(self, parent):
        """Create pipeline status indicator"""
        status_frame = tk.Frame(parent, bg=ModernStyle.COLORS['bg_light'])
        status_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.status_var = tk.StringVar(value="Ready to configure and run pipeline")
        status_label = tk.Label(status_frame, textvariable=self.status_var,
                              font=ModernStyle.get_font('body'),
                              fg=ModernStyle.COLORS['success'],
                              bg=ModernStyle.COLORS['bg_light'])
        status_label.pack(anchor='w')

    def create_configuration_tab(self):
        """Enhanced configuration tab"""
        self.config_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.config_tab, text="Configuration")

        # Scrollable frame for better organization
        canvas = tk.Canvas(self.config_tab, bg=ModernStyle.COLORS['bg_light'])
        scrollbar = ttk.Scrollbar(self.config_tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Create parameter widgets
        self.create_parameter_widgets(scrollable_frame)

    def create_advanced_tab(self):
        """Enhanced advanced parameters tab with crisp text"""
        # Use ttk.Frame for the tab container
        self.advanced_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.advanced_tab, text="Advanced")

        # Outer frame inside tab (for padding)
        outer_frame = ttk.Frame(self.advanced_tab, style="Card.TFrame")
        outer_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)

        # Create advanced widgets
        self.create_advanced_widgets(outer_frame)

        # Steps overview section using ttk.LabelFrame
        steps_frame = ttk.LabelFrame(
            outer_frame,
            text="Pipeline Steps Overview",
            style="Card.TFrame"  # Background color comes from Card.TFrame
        )
        steps_frame.pack(fill=tk.BOTH, expand=True, pady=20)

        # Add ScrolledText (tk.Text is necessary, but font/bg can be scaled)
        steps_text = scrolledtext.ScrolledText(
            steps_frame,
            height=15,
            font=ModernStyle.get_font('code'),  # Make sure 'code' font is DPI-scaled
            background=ModernStyle.COLORS['bg_card'],
            foreground="#222222",
            borderwidth=0,
            highlightthickness=0
        )
        steps_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        steps_text.insert(tk.END, self._get_pipeline_steps_text())
        steps_text.config(state=tk.DISABLED)

    def create_monitoring_tab(self):
        """Enhanced monitoring tab"""
        self.monitor_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.monitor_tab, text="Monitoring")

        frame = ttk.Frame(self.monitor_tab, style="Card.TFrame")
        frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)

        # Progress section
        progress_frame = ttk.LabelFrame(frame, text="Pipeline Progress", 
                               style="Large.TLabelframe")
        progress_frame.pack(fill=tk.X, pady=(0, 20))
        

        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(progress_frame, variable=self.progress_var, 
                                           mode='indeterminate', length=800)
        self.progress_bar.pack(pady=30, padx=20)

        # Log controls
        log_controls = ttk.Frame(frame, style="Card.TFrame")
        log_controls.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(log_controls, text="Clear Log",
          command=self.clear_log, style="Primary.TButton").pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(log_controls, text="Save Log",
          command=self.save_log, style="Primary.TButton").pack(side=tk.LEFT)

        # Enhanced log window
        log_frame = ttk.LabelFrame(frame, text="Pipeline Output Log",
                          style="Card.TFrame",)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = scrolledtext.ScrolledText(log_frame, 
                                                 font=ModernStyle.get_font('code'),
                                                 bg=ModernStyle.COLORS['bg_dark'], 
                                                 fg='#00FF7F',
                                                 insertbackground='#00FF7F')
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Initial log messages
        self.log_message("MSGFPDOWN Pipeline GUI Enhanced - Ready")
        if self.pipeline_path:
            self.log_message(f"Pipeline detected: {self.pipeline_path}")
    
    def create_help_tab(self):
        """Enhanced help tab"""
        self.help_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.help_tab, text="Help")

        frame = tk.Frame(self.help_tab, bg=ModernStyle.COLORS['bg_light'])
        frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)

        help_text = scrolledtext.ScrolledText(frame, font=ModernStyle.get_font('body'),
                                             bg=ModernStyle.COLORS['bg_card'])
        help_text.pack(fill=tk.BOTH, expand=True)
        help_text.insert(tk.END, self._get_help_text())
        help_text.config(state=tk.DISABLED)

    def create_footer(self, parent):
        """Enhanced footer with better button layout"""
        footer = tk.Frame(parent, bg="#f0f0f0")
        footer.pack(fill=tk.X, pady=(20, 0))

        # Left controls - main actions
        left_controls = tk.Frame(footer, bg="#f0f0f0")
        left_controls.pack(side=tk.LEFT)

        self.run_button = ttk.Button(left_controls, text="Run Pipeline", 
                                    command=self.run_pipeline, style='Success.TButton')
        self.run_button.pack(side=tk.LEFT, padx=5)
        EnhancedToolTip(self.run_button, "Start the MSGFPDOWN pipeline analysis")

        self.stop_button = ttk.Button(left_controls, text="Stop Pipeline", 
                                     command=self.stop_pipeline, style='Danger.TButton', 
                                     state=tk.DISABLED)
        self.stop_button.pack(side=tk.LEFT, padx=5)
        EnhancedToolTip(self.stop_button, "Terminate the running pipeline")

        # Right controls - utility actions
        right_controls = tk.Frame(footer, bg=ModernStyle.COLORS['bg_light'])
        right_controls.pack(side=tk.RIGHT)

        version_text = "v1.3"
        tk.Label(
            right_controls,
            text=version_text,
            font=ModernStyle.get_font('small'),
            bg=ModernStyle.COLORS['bg_light']
        ).pack(side=tk.RIGHT, padx=(10, 0))

        ttk.Button(right_controls, text="Quit", 
                  command=self.on_close).pack(side=tk.LEFT)

    # -------------------- Helper Methods for UI Creation --------------------
    def create_section(self, parent, title):
        section_frame = ttk.Frame(parent, style="Card.TFrame")
        section_frame.pack(fill=tk.X, pady=10, padx=10)

        # Use ttk.Label for title
        title_label = ttk.Label(section_frame, text=title, style="Subtitle.TLabel")
        title_label.pack(anchor=tk.W, pady=(0,5))

        return section_frame

    def create_parameter_widgets(self, parent):
        """Create enhanced parameter input widgets"""
        # Required parameters section
        req_section = self.create_section(parent, "Required Parameters")
        
        self.create_file_input(req_section, "SIC Input Directory:", 
                              self.sic_input_dir, self.browse_sic_input_dir)
        self.create_file_input(req_section, "Final Output File:", 
                              self.final_output, self.browse_final_output)  
        self.create_file_input(req_section, "Database Directory:", 
                              self.database_path, self.browse_database_dir)

        # Filtering parameters section
        filter_section = self.create_section(parent, "Filtering Parameters")
        
        # Score fields in a grid
        grid_frame = tk.Frame(filter_section, bg=ModernStyle.COLORS['bg_card'])
        grid_frame.pack(fill=tk.X, padx=10, pady=10)

        # Primary score
        tk.Label(grid_frame, text="Primary Score Field:", 
                font=ModernStyle.get_font('body')).grid(row=0, column=0, sticky='w', padx=(0,10))
        ttk.Entry(grid_frame, textvariable=self.score_field, 
                 width=20).grid(row=0, column=1, padx=(0,20))
        tk.Label(grid_frame, text="Threshold:", 
                font=ModernStyle.get_font('body')).grid(row=0, column=2, sticky='w', padx=(0,10))
        ttk.Entry(grid_frame, textvariable=self.threshold, 
                 width=10).grid(row=0, column=3)

        # Secondary score
        tk.Label(grid_frame, text="Secondary Score Field:", 
                font=ModernStyle.get_font('body')).grid(row=1, column=0, sticky='w', padx=(0,10), pady=(10,0))
        ttk.Entry(grid_frame, textvariable=self.score_field2, 
                 width=20).grid(row=1, column=1, padx=(0,20), pady=(10,0))
        tk.Label(grid_frame, text="Threshold:", 
                font=ModernStyle.get_font('body')).grid(row=1, column=2, sticky='w', padx=(0,10), pady=(10,0))
        ttk.Entry(grid_frame, textvariable=self.threshold2, 
                 width=10).grid(row=1, column=3, pady=(10,0))

        # Analysis parameters section
        analysis_section = self.create_section(parent, "Analysis Parameters")
        
        analysis_grid = tk.Frame(analysis_section, bg=ModernStyle.COLORS['bg_card'])
        analysis_grid.pack(fill=tk.X, padx=10, pady=10)

        # Min peptides
        tk.Label(analysis_grid, text="Min Peptides:", 
                font=ModernStyle.get_font('body')).grid(row=0, column=0, sticky='w', padx=(0,10))
        ttk.Spinbox(analysis_grid, from_=1, to=50, textvariable=self.num_pep, 
                   width=8).grid(row=0, column=1, padx=(0,20))

        # Mode selection
        tk.Label(analysis_grid, text="Mode:", 
                font=ModernStyle.get_font('body')).grid(row=0, column=2, sticky='w', padx=(0,10))
        mode_combo = ttk.Combobox(analysis_grid, textvariable=self.mode,
                                 values=["unique_only", "requires_unique", "all_matches"],
                                 width=18, state='readonly')
        mode_combo.grid(row=0, column=3)

        # Rollup method
        tk.Label(analysis_grid, text="Rollup Method:", 
                font=ModernStyle.get_font('body')).grid(row=1, column=0, sticky='w', padx=(0,10), pady=(10,0))
        rollup_combo = ttk.Combobox(analysis_grid, textvariable=self.rollup,
                                   values=["sum", "rrollup", "zrollup", "qrollup", "weighted"],
                                   width=18, state='readonly')
        rollup_combo.grid(row=1, column=1, columnspan=2, sticky='w', pady=(10,0))

        # Cluster step checkbox
        cluster_frame = tk.Frame(analysis_section, bg=ModernStyle.COLORS['bg_card'])
        cluster_frame.pack(fill=tk.X, padx=10, pady=10)
        ttk.Checkbutton(cluster_frame, text="Run Cluster Crosstab Generation (Step 10)", 
                       variable=self.run_cluster_step).pack(side=tk.LEFT)

    def create_advanced_widgets(self, parent):
        """Create advanced parameter widgets"""
        advanced_section = self.create_section(parent, "Advanced Parameters")
        
        # Outlier alpha
        outlier_frame = tk.Frame(advanced_section, bg=ModernStyle.COLORS['bg_card'])
        outlier_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Label(outlier_frame, text="Outlier Alpha (RRrollup only):", 
                font=ModernStyle.get_font('body')).pack(side=tk.LEFT)
        ttk.Entry(outlier_frame, textvariable=self.outlier_alpha, 
                 width=12).pack(side=tk.LEFT, padx=(10,0))

        # Coverage TSV
        self.create_file_input(advanced_section, "Coverage TSV (optional):", 
                              self.coverage_tsv, self.browse_coverage_tsv)
        
        # Cluster Crosstab TSV - New addition
        self.create_file_input(advanced_section, "Cluster Output (Required when Cluster Rollup opted):", 
                              self.cluster_crosstab_tsv, self.browse_cluster_crosstab_tsv)

    def create_file_input(self, parent, label_text, variable, browse_command, tooltip=""):
        """Create file input widget with browse button"""
        input_frame = tk.Frame(parent, bg=ModernStyle.COLORS['bg_card'])
        input_frame.pack(fill=tk.X, pady=8, padx=10)
        
        tk.Label(input_frame, text=label_text, 
                font=ModernStyle.get_font('body')).pack(anchor='w')
        
        row_frame = tk.Frame(input_frame, bg=ModernStyle.COLORS['bg_card'])
        row_frame.pack(fill=tk.X, pady=(5,0))
        
        entry = ttk.Entry(row_frame, textvariable=variable, style='Modern.TEntry')
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0,10))
        
        btn = ttk.Button(row_frame, text="Browse", command=browse_command)
        btn.pack(side=tk.RIGHT)
        
        if tooltip:
            EnhancedToolTip(btn, tooltip)

    # -------------------- Browse Methods --------------------
    def browse_sic_input_dir(self):
        """Browse for SIC input directory"""
        directory = filedialog.askdirectory(title="Select SIC Input Directory")
        if directory:
            self.sic_input_dir.set(directory)

    def browse_final_output(self):
        """Browse for final output file"""
        filename = filedialog.asksaveasfilename(
            title="Choose Final Output File",
            defaultextension=".tsv",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if filename:
            self.final_output.set(filename)

    def browse_database_dir(self):
        """Browse for database directory"""
        directory = filedialog.askdirectory(title="Select Database Directory")
        if directory:
            self.database_path.set(directory)

    def browse_coverage_tsv(self):
        """Browse for coverage TSV file"""
        filename = filedialog.askopenfilename(
            title="Select Coverage TSV File",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if filename:
            self.coverage_tsv.set(filename)
    
    #clsuter output - New addition
    def browse_cluster_crosstab_tsv(self):
        """Browse for Cluster Crosstab TSV file"""
        filename = filedialog.asksaveasfilename(
            title="Select Cluster Crosstab TSV File",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")]
        )
        if filename:
            self.cluster_crosstab_tsv.set(filename)

    # -------------------- Pipeline Execution --------------------
    def run_pipeline(self):
        """Execute the pipeline with enhanced error handling"""
        if self.is_running:
            messagebox.showinfo("Pipeline Running", "Pipeline is already running.")
            return

        # Validate required inputs
        if not self.sic_input_dir.get() or not Path(self.sic_input_dir.get()).is_dir():
            messagebox.showerror("Missing Input", "Please select a valid SIC input directory.")
            self.notebook.select(0)  # Switch to config tab
            return

        if not self.final_output.get():
            messagebox.showerror("Missing Output", "Please specify a final output file path.")
            self.notebook.select(0)
            return

        if not self.database_path.get() or not Path(self.database_path.get()).is_dir():
            messagebox.showerror("Missing Database", "Please select a valid database directory.")
            self.notebook.select(0)
            return

        # Prepare command arguments
        # if self.pipeline_path.startswith("python -m "):
        #     # Handle module execution
        #     module_name = self.pipeline_path.split("python -m ")[1]
        #     cmd = [sys.executable, "-m", module_name]
        # else:
        #     # Handle direct script execution
        #     cmd = [sys.executable, self.pipeline_path]
        cmd = [sys.executable, self.pipeline_path]
        cmd.extend(["-i", self.sic_input_dir.get()])
        cmd.extend(["-o", self.final_output.get()])
        cmd.extend(["--database", self.database_path.get()])

        # Add optional parameters
        if self.score_field.get():
            cmd.extend(["--score_field", self.score_field.get()])
        if self.threshold.get():
            cmd.extend(["--threshold", self.threshold.get()])
        if self.score_field2.get():
            cmd.extend(["--score_field2", self.score_field2.get()])
        if self.threshold2.get():
            cmd.extend(["--threshold2", self.threshold2.get()])
        if self.num_pep.get():
            cmd.extend(["--num_pep", self.num_pep.get()])
        if self.mode.get():
            cmd.extend(["--mode", self.mode.get()])
        if self.rollup.get():
            cmd.extend(["--rollup", self.rollup.get()])
        if self.outlier_alpha.get():
            cmd.extend(["--outlier_alpha", self.outlier_alpha.get()])
        if self.coverage_tsv.get():
            cmd.extend(["--coverage_tsv", self.coverage_tsv.get()])
        if self.cluster_crosstab_tsv.get():
            cmd.extend(["--cluster_output", self.cluster_crosstab_tsv.get()]) #New addition
        if self.run_cluster_step.get():
            cmd.append("--cluster_rollup")

        # Start pipeline execution
        self.log_message("="*50)
        self.log_message("Starting MSGFPDOWN Pipeline")
        self.log_message(f"Command: {' '.join(cmd[:6])} ...")
        self.log_message("="*50)

        self.status_var.set("Pipeline starting...")
        self.is_running = True
        self._stop_requested = False
        self.run_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)

        # Switch to monitoring tab
        self.notebook.select(2)

        try:
            self.process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                bufsize=1, 
                text=True, 
                universal_newlines=True
            )
        except Exception as e:
            self.log_message(f"Failed to start pipeline: {e}")
            self.status_var.set("Failed to start")
            self.is_running = False
            self.run_button.config(state=tk.NORMAL)
            self.stop_button.config(state=tk.DISABLED)
            return

        # Start progress animation
        self.progress_bar.start(10)

        # Start monitoring threads
        threading.Thread(target=self._reader_thread, daemon=True).start()
        threading.Thread(target=self._monitor_thread, daemon=True).start()

    def stop_pipeline(self):
        """Stop the running pipeline"""
        if not self.is_running or not self.process:
            return

        if not messagebox.askyesno("Stop Pipeline", 
                                  "Are you sure you want to stop the running pipeline?\n"
                                  "This will terminate the current analysis."):
            return

        self._stop_requested = True
        self.status_var.set("Stopping pipeline...")
        self.log_message("Stop requested - terminating pipeline...")

        try:
            self.process.terminate()
        except Exception as e:
            self.log_message(f"Error during termination: {e}")

        # Start timeout watcher
        threading.Thread(target=self._terminate_watcher, daemon=True).start()

    def _reader_thread(self):
        """Thread to read pipeline output"""
        try:
            if self.process and self.process.stdout:
                for line in iter(self.process.stdout.readline, ''):
                    if line:
                        self._log_queue.put(line.rstrip())
                self.process.stdout.close()
        except Exception as e:
            self._log_queue.put(f"[Reader Error] {e}")

    def _monitor_thread(self):
        """Thread to monitor pipeline completion"""
        try:
            return_code = self.process.wait()
            if return_code == 0:
                self._log_queue.put("[SUCCESS] Pipeline completed successfully!")
            else:
                self._log_queue.put(f"[ERROR] Pipeline failed with exit code {return_code}")
        except Exception as e:
            self._log_queue.put(f"[Monitor Error] {e}")
        finally:
            self.is_running = False
            self.root.after(0, self._on_process_end)

    def _on_process_end(self):
        """Handle pipeline completion"""
        self.progress_bar.stop()
        self.run_button.config(state=tk.NORMAL)
        self.stop_button.config(state=tk.DISABLED)
        
        if self._stop_requested:
            self.status_var.set("Pipeline stopped by user")
        else:
            self.status_var.set("Pipeline completed")

    def _terminate_watcher(self, timeout=10):
        """Watch for process termination and force kill if needed"""
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.process.poll() is not None:
                return
            time.sleep(0.2)

        # Force kill if still running
        try:
            self.log_message("Process unresponsive - force killing...")
            self.process.kill()
        except Exception as e:
            self.log_message(f"Failed to force kill: {e}")

    # -------------------- Logging Methods --------------------
    def log_message(self, message):
        """Add message to log with timestamp"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def _poll_log_queue(self):
        """Poll log queue for new messages"""
        try:
            while not self._log_queue.empty():
                message = self._log_queue.get_nowait()
                timestamp = datetime.now().strftime("%H:%M:%S")
                self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
                self.log_text.see(tk.END)
        except queue.Empty:
            pass
        except Exception as e:
            print(f"Log queue error: {e}")
        finally:
            self.root.after(200, self._poll_log_queue)

    def clear_log(self):
        """Clear the log window"""
        if messagebox.askyesno("Clear Log", "Clear the current log output?"):
            self.log_text.delete("1.0", tk.END)
            self.log_message("Log cleared")

    def save_log(self):
        """Save log to file"""
        filename = filedialog.asksaveasfilename(
            title="Save Log As...",
            defaultextension=".log",
            filetypes=[("Log files", "*.log"), ("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            try:
                with open(filename, "w", encoding="utf-8") as f:
                    f.write(self.log_text.get("1.0", tk.END))
                messagebox.showinfo("Log Saved", f"Log saved to {filename}")
            except Exception as e:
                messagebox.showerror("Save Error", f"Failed to save log: {e}")

    # -------------------- Help and Information --------------------
    def _get_pipeline_steps_text(self):
        """Get pipeline steps description"""
        return """MSGFPDOWN Pipeline Steps Overview

Step 1: Prepare SICs
    - Processes input SIC files for analysis
    - Input: Raw SIC directory
    - Output: Processed SIC files

Step 2: FDR Estimation  
    - Performs False Discovery Rate estimation
    - Input: Processed SIC files
    - Output: FDR estimated files

Step 3: Merge Shared PHRP
    - Merges shared PHRP data
    - Input: FDR estimated files
    - Output: Merged PHRP data

Step 4: True PEP Filter
    - Filters peptides based on PEP scores and thresholds
    - Uses: score_field, threshold, filter_method parameters
    - Input: Merged PHRP data
    - Output: Filtered peptide data

Step 5: Peptide Crosstab Generation
    - Creates peptide crosstab matrix
    - Input: Filtered peptide data
    - Output: Peptide crosstab TSV

Step 6: FASTA Cluster Generator
    - Generates protein sequence clusters
    - Input: Database files
    - Output: Clustered sequences

Step 7: Peptide Crosstab Annotator
    - Annotates peptide crosstab with protein information
    - Uses: database_path parameter
    - Input: Peptide crosstab + Database
    - Output: Annotated peptide crosstab

Step 8: Peptide-Protein Map Builder
    - Builds peptide to protein mapping
    - Uses: num_pep, mode parameters
    - Input: Annotated peptide crosstab
    - Output: Peptide-protein map

Step 9: Protein Rollup
    - Rolls peptide data up to protein level
    - Uses: rollup, mode, outlier_alpha, coverage_tsv parameters
    - Input: Peptide-protein map
    - Output: Final protein results

Step 10: Cluster Crosstab Generator (Optional)
    - Generates cluster-level crosstab
    - Enabled by: run_cluster_step parameter
    - Input: Annotated data
    - Output: Cluster crosstab
"""

    def _get_help_text(self):
        """Get help text for the application"""
        return """MSGFPDOWN Pipeline Enhanced GUI - Help

OVERVIEW:
This GUI provides an enhanced interface for running the MSGFPDOWN proteomics pipeline.
The pipeline automatically detects the script location and provides a modern,
high-DPI compatible interface.

REQUIRED PARAMETERS:
- SIC Input Directory: Directory containing your SIC files for analysis
- Final Output File: Where you want the final results saved (.tsv format)
- Database Directory: Directory containing protein database files

OPTIONAL PARAMETERS:
- Primary/Secondary Score Fields: Score fields for peptide filtering (e.g., MSMSScore)
- Thresholds: Numeric thresholds for score-based filtering
- Filter Method: Specific filtering method for Step 4
- Min Peptides: Minimum peptides required per protein (default: 2)
- Mode: Protein inference mode (unique_only, requires_unique, all_matches)
- Rollup Method: How to roll peptides to proteins (sum, rrollup, zrollup, qrollup, weighted)

ADVANCED PARAMETERS:
- Outlier Alpha: Alpha value for Grubbs outlier detection (RRrollup only)
- Coverage TSV: Optional coverage file for analysis
- Run Cluster Step: Enable optional Step 10 (cluster crosstab generation)

USAGE TIPS:
1. Fill in required parameters first
2. Configure optional parameters as needed
3. Use the Monitoring tab to watch pipeline progress
4. Logs are automatically captured and can be saved
5. Pipeline can be stopped gracefully if needed

KEYBOARD SHORTCUTS:
- Ctrl+R: Run pipeline
- Ctrl+S: Stop pipeline
- Ctrl+Q: Quit application
- F1: Show help tab

SYSTEM REQUIREMENTS:
- Python environment with MSGFPDOWN pipeline installed
- Cross-platform compatible (Linux, macOS, Windows)
- High-DPI display support for crisp rendering

For technical support, refer to the MSGFPDOWN pipeline documentation.
"""

    # -------------------- Keyboard Shortcuts --------------------
    def bind_shortcuts(self):
        """Bind keyboard shortcuts"""
        self.root.bind("<Control-r>", lambda e: self.run_pipeline() if not self.is_running else None)
        self.root.bind("<Control-s>", lambda e: self.stop_pipeline() if self.is_running else None)
        self.root.bind("<Control-q>", lambda e: self.on_close())
        self.root.bind("<F1>", lambda e: self.notebook.select(3))  # Help tab

    def on_close(self):
        """Handle application close"""
        if self.is_running:
            if messagebox.askyesno("Quit Application", 
                                  "Pipeline is running. Quit anyway?\n"
                                  "This will terminate the pipeline."):
                try:
                    if self.process:
                        self.process.terminate()
                        time.sleep(0.5)
                        if self.process.poll() is None:
                            self.process.kill()
                except Exception:
                    pass
            else:
                return

        self.root.destroy()

# -------------------- Application Entry Point --------------------
def main():
    print("Launching MS-GF+ Downstream GUI")
    """Main application entry point"""
    root = tk.Tk()
    app = MSGFPDOWNEnhancedGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()