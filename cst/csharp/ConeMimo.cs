// Конус MIMO: шаг и полоса. Изображение точечной цели по методу Чапурского для конусной MIMO-системы
// (4 кольца: приёмники на 1-м и 3-м, передатчики на 2-м и 4-м), цель далеко на горизонте в азимуте 0°.
// Расчёт: Ψ(u) = |Σ_f AF_tx(f,u)·AF_rx(f,u)|, AF = Σ exp[jk r·(u − u0)], частоты через 0,25 ГГц вокруг 10 ГГц.
// Сборка (.NET Framework 4, есть в любой Windows 10/11):
//   C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe /codepage:65001 /target:winexe /optimize+ /win32icon:radar.ico /out:ConeMimo.exe ConeMimo.cs
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Globalization;
using System.Linq;
using System.Windows.Forms;

static class Calc
{
    public const double C = 299.792458;             // мм·ГГц
    public const double LAM = C / 10.0;             // λ на 10 ГГц, мм
    public static readonly double[][] RINGS = {     // R/λ, z/λ, n, 1 = передатчики
        new double[] { 5.0, 0.0, 10, 0 }, new double[] { 4.17, 1.4, 8, 1 },
        new double[] { 3.33, 2.8, 8, 0 }, new double[] { 2.5, 4.2, 6, 1 } };
    public static readonly double D0 = 2 * Math.PI * 5.0 / 10;   // шаг основного конуса, 3,14λ

    public class Geo { public List<double[]> Rx = new List<double[]>(), Tx = new List<double[]>(); public double S, R1, RingStep; }

    public static Geo Geometry(double d, bool scale)
    {
        var g = new Geo(); g.S = scale ? d / D0 : 1.0;
        foreach (var r in RINGS)
        {
            double Rm = r[0] * LAM * g.S, zm = r[1] * LAM * g.S;
            int n = scale ? (int)r[2] : Math.Max(3, (int)Math.Round(2 * Math.PI * Rm / (d * LAM)));
            for (int i = 0; i < n; i++)
            {
                double a = 2 * Math.PI * i / n + (r[3] > 0 ? Math.PI / n : 0);
                (r[3] > 0 ? g.Tx : g.Rx).Add(new[] { Rm * Math.Cos(a), Rm * Math.Sin(a), zm });
            }
        }
        g.R1 = 5 * LAM * g.S; g.RingStep = Math.Sqrt(0.83 * 0.83 + 1.4 * 1.4) * LAM * g.S;
        return g;
    }

    public static double[] Freqs(int bPct)
    {
        double half = 10.0 * bPct / 200; int n = (int)Math.Round(2 * half / 0.25);
        var f = new double[n + 1];
        for (int i = 0; i <= n; i++) f[i] = 10 - half + 0.25 * i;
        return f;
    }

    static void Af(List<double[]> pts, double k, double ux, double uy, double uz, out double re, out double im)
    {
        re = 0; im = 0;
        foreach (var p in pts)
        {
            double ph = k * (p[0] * (ux - 1) + p[1] * uy + p[2] * uz);
            re += Math.Cos(ph); im += Math.Sin(ph);
        }
    }

    public static double[] Cut(Geo g, double[] fs, double[][] dirs)
    {
        var o = new double[dirs.Length]; double mx = 0;
        for (int i = 0; i < dirs.Length; i++)
        {
            double re = 0, im = 0;
            foreach (double f in fs)
            {
                double k = 2 * Math.PI * f / C, a, b, c, d;
                Af(g.Tx, k, dirs[i][0], dirs[i][1], dirs[i][2], out a, out b);
                Af(g.Rx, k, dirs[i][0], dirs[i][1], dirs[i][2], out c, out d);
                re += a * c - b * d; im += a * d + b * c;
            }
            o[i] = Math.Sqrt(re * re + im * im); if (o[i] > mx) mx = o[i];
        }
        for (int i = 0; i < o.Length; i++) o[i] = 20 * Math.Log10(o[i] / mx + 1e-12);
        return o;
    }

    public static double[] Range(double[] fs, double[] dr)
    {
        var o = new double[dr.Length]; double mx = 0;
        for (int i = 0; i < dr.Length; i++)
        {
            double re = 0, im = 0;
            foreach (double f in fs) { double ph = 4 * Math.PI * f / C * dr[i]; re += Math.Cos(ph); im += Math.Sin(ph); }
            o[i] = Math.Sqrt(re * re + im * im); if (o[i] > mx) mx = o[i];
        }
        for (int i = 0; i < o.Length; i++) o[i] = 20 * Math.Log10(o[i] / mx + 1e-12);
        return o;
    }

    // ширина главного лепестка по −3 дБ и наибольший уровень вне главного лепестка (до первых минимумов)
    public static void Lobe(double[] x, double[] d, bool circ, out double bw, out double sll)
    {
        int n = d.Length, i = 0;
        for (int j = 1; j < n; j++) if (d[j] > d[i]) i = j;
        Func<int, double> at = j => circ ? d[((j % n) + n) % n] : d[Math.Min(Math.Max(j, 0), n - 1)];
        int lo = 0, hi = 0;
        while (lo > -n / 2 && at(i + lo - 1) >= -3 && (circ || i + lo > 0)) lo--;
        while (hi < n / 2 && at(i + hi + 1) >= -3 && (circ || i + hi < n - 1)) hi++;
        int l = 0, r = 0;
        while (l > -n / 2 && at(i + l - 1) < at(i + l) && (circ || i + l > 0)) l--;
        while (r < n / 2 && at(i + r + 1) < at(i + r) && (circ || i + r < n - 1)) r++;
        var main = new HashSet<int>();
        for (int j = l; j <= r; j++) main.Add(((i + j) % n + n) % n);
        sll = double.NegativeInfinity;
        for (int j = 0; j < n; j++) if (!main.Contains(j) && d[j] > sll) sll = d[j];
        bool flat = !circ && i + lo <= 0 && i + hi >= n - 1;
        bw = flat ? double.NaN : (hi - lo) * (x[1] - x[0]);
    }
}

class Plot : Panel
{
    public string Title = "", XLabel = "", YLabel = "Ψ, дБ";
    public double X0, X1, Y0 = -40, Y1 = 0; public double[] XT, YT = { -40, -30, -20, -10, 0 };
    public double[] X; public List<Tuple<double[], Color, bool>> Series = new List<Tuple<double[], Color, bool>>();
    public List<Tuple<double, Color>> VLines = new List<Tuple<double, Color>>(); public string VLabel;
    public List<Tuple<string, Color, bool>> Legend = new List<Tuple<string, Color, bool>>();   // подпись, цвет, пунктир
    public Plot() { DoubleBuffered = true; BackColor = Color.White; ResizeRedraw = true; }

    protected override void OnPaint(PaintEventArgs e)
    {
        var g = e.Graphics; g.SmoothingMode = SmoothingMode.AntiAlias;
        float q = g.DpiX / 96f; int L = (int)(50 * q), R = (int)(12 * q), T = (int)(26 * q), B = (int)(38 * q); float pw = Width - L - R, ph = Height - T - B;
        if (pw < 20 || ph < 20) return;
        Func<double, float> PX = v => (float)(L + (v - X0) / (X1 - X0) * pw);
        Func<double, float> PY = v => (float)(T + (Y1 - Math.Max(Y0, Math.Min(Y1, v))) / (Y1 - Y0) * ph);
        using (var fTitle = new Font("Segoe UI", 9.5f, FontStyle.Bold))
        using (var f = new Font("Segoe UI", 8.5f))
        using (var gridPen = new Pen(Color.FromArgb(228, 232, 237)))
        using (var fr = new Pen(Color.FromArgb(190, 198, 208)))
        {
            g.DrawString(Title, fTitle, Brushes.Black, 6, 4);
            // легенда в строке заголовка, прижата вправо: образец линии и подпись
            {
                float x = Width - 8 * q, yl = 4 * q + g.MeasureString("Ag", fTitle).Height / 2;
                for (int i = Legend.Count - 1; i >= 0; i--)
                {
                    var it = Legend[i]; var sz = g.MeasureString(it.Item1, f);
                    x -= sz.Width; g.DrawString(it.Item1, f, Brushes.DimGray, x, yl - sz.Height / 2);
                    x -= 24 * q;
                    using (var lp = new Pen(it.Item2, it.Item3 ? 1.3f : 2f)) { if (it.Item3) lp.DashPattern = new float[] { 4, 3 }; g.DrawLine(lp, x, yl, x + 20 * q, yl); }
                    x -= 14 * q;
                }
            }
            var sc = new StringFormat { Alignment = StringAlignment.Center };
            var sr = new StringFormat { Alignment = StringAlignment.Far, LineAlignment = StringAlignment.Center };
            foreach (double t in XT) { g.DrawLine(gridPen, PX(t), T, PX(t), T + ph); g.DrawString(Num(t), f, Brushes.DimGray, PX(t), T + ph + 3 * q, sc); }
            foreach (double t in YT) { g.DrawLine(gridPen, L, PY(t), L + pw, PY(t)); g.DrawString(Num(t), f, Brushes.DimGray, L - 4, PY(t), sr); }
            g.DrawString(XLabel, f, Brushes.DimGray, L + pw / 2, T + ph + 19 * q, sc);
            var st = g.Save(); g.TranslateTransform(12 * q, T + ph / 2); g.RotateTransform(-90); g.DrawString(YLabel, f, Brushes.DimGray, 0, -8, sc); g.Restore(st);
            g.DrawRectangle(fr, L, T, pw, ph);
            g.SetClip(new RectangleF(L, T, pw, ph));
            foreach (var v in VLines) using (var p = new Pen(v.Item2, 1.2f) { DashStyle = DashStyle.Dash }) g.DrawLine(p, PX(v.Item1), T, PX(v.Item1), T + ph);
            if (VLines.Count > 0 && VLabel != null)
                using (var b = new SolidBrush(VLines[0].Item2))
                    g.DrawString(VLabel, f, b, PX(VLines[VLines.Count - 1].Item1) + 3 * q, T + 2 * q);
            // уровень минус 3 дБ (по нему считаются ширина луча и разрешение по дальности): тонкий пунктир цвета кривой
            using (var p = new Pen(Color.FromArgb(150, 15, 107, 92), 1f) { DashPattern = new float[] { 5, 4 } })
                g.DrawLine(p, L, PY(-3), L + pw, PY(-3));
            foreach (var s in Series)
            {
                if (s.Item1 == null || X == null) continue;
                var pts = new PointF[X.Length];
                for (int i = 0; i < X.Length; i++) pts[i] = new PointF(PX(X[i]), PY(s.Item1[i]));
                using (var p = new Pen(s.Item2, s.Item3 ? 1.2f : 1.7f)) { if (s.Item3) p.DashStyle = DashStyle.Dash; g.DrawLines(p, pts); }
            }
            g.ResetClip();
            // подпись уровня на оси: короткая засечка и «−3» в цвет пунктира
            using (var ac = new SolidBrush(Color.FromArgb(15, 107, 92)))
            using (var fb = new Font("Segoe UI", 8.5f, FontStyle.Bold))
            using (var tp = new Pen(Color.FromArgb(15, 107, 92), 1.5f))
            {
                float y3 = PY(-3); g.DrawLine(tp, L - 4 * q, y3, L, y3);
                g.DrawString("−3", fb, ac, L - 5 * q, y3, sr);
            }
        }
    }
    static string Num(double v) { return v.ToString("0.##", CultureInfo.InvariantCulture).Replace("-", "−").Replace(".", ","); }
}

class GeoView : Panel
{
    public Calc.Geo G;
    public GeoView() { DoubleBuffered = true; BackColor = Color.White; ResizeRedraw = true; }
    protected override void OnPaint(PaintEventArgs e)
    {
        var g = e.Graphics; g.SmoothingMode = SmoothingMode.AntiAlias;
        using (var fTitle = new Font("Segoe UI", 9.5f, FontStyle.Bold)) g.DrawString("Расстановка, вид сверху", fTitle, Brushes.Black, 6, 4);
        if (G == null) return;
        double Rm = 0; foreach (var p in G.Rx.Concat(G.Tx)) Rm = Math.Max(Rm, Math.Sqrt(p[0] * p[0] + p[1] * p[1]));
        // постоянный масштаб в миллиметрах: рамка рассчитана на самый большой конус (нижнее кольцо 150 мм),
        // поэтому при малом шаге в режиме масштаба конус на картинке действительно маленький
        const double RMAX = 5.0 * Calc.LAM;
        float q0 = e.Graphics.DpiX / 96f;
        float cx = Width / 2f, cy = (Height + 20 * q0 - 24 * q0) / 2f + 6 * q0;
        float sc = (float)((Math.Min(Width, Height - 60 * q0) / 2.0 - 10 * q0) / RMAX);
        using (var ring = new Pen(Color.FromArgb(228, 232, 237)))
            foreach (var r in Calc.RINGS) { float rr = (float)(r[0] * Calc.LAM * G.S * sc); g.DrawEllipse(ring, cx - rr, cy - rr, 2 * rr, 2 * rr); }
        using (var vb = new SolidBrush(Color.FromArgb(150, 160, 172)))
        {
            int step = Math.Max(1, G.Rx.Count * G.Tx.Count / 3000), c = 0;
            foreach (var t in G.Tx) foreach (var r in G.Rx)
                { if (c++ % step != 0) continue; g.FillRectangle(vb, cx + (float)((t[0] + r[0]) / 2 * sc) - 0.8f, cy - (float)((t[1] + r[1]) / 2 * sc) - 0.8f, 1.6f, 1.6f); }
        }
        using (var bRx = new SolidBrush(Color.FromArgb(29, 95, 168)))
        using (var bTx = new SolidBrush(Color.FromArgb(194, 65, 12)))
        {
            foreach (var p in G.Rx) g.FillEllipse(bRx, cx + (float)(p[0] * sc) - 3, cy - (float)(p[1] * sc) - 3, 6, 6);
            foreach (var p in G.Tx) g.FillEllipse(bTx, cx + (float)(p[0] * sc) - 3, cy - (float)(p[1] * sc) - 3, 6, 6);
            using (var f = new Font("Segoe UI", 8.5f))
            {
                float q = g.DpiX / 96f, y = Height - 20 * q, x = 8 * q;
                Action<Brush, bool, string> item = (br, round, txt) =>
                {
                    float m = 7 * q, yy = y + 4 * q;
                    if (round) g.FillEllipse(br, x, yy, m, m); else g.FillRectangle(br, x + 2 * q, yy + 2 * q, 3 * q, 3 * q);
                    x += m + 4 * q; g.DrawString(txt, f, Brushes.DimGray, x, y);
                    x += g.MeasureString(txt, f).Width + 12 * q;
                };
                item(bRx, true, "приёмники"); item(bTx, true, "передатчики"); item(Brushes.Gray, false, "виртуальные центры");
            }
        }
        using (var p = new Pen(Color.FromArgb(15, 107, 92), 1.5f)) g.DrawLine(p, cx + (float)(Rm * sc) + 4, cy, Width - 6, cy);
        // масштабная линейка 50 мм в правом нижнем углу
        using (var p = new Pen(Color.FromArgb(60, 70, 80), 1.5f))
        using (var f = new Font("Segoe UI", 8.5f))
        {
            float L = (float)(50 * sc), x1 = Width - 10 * q0, x0 = x1 - L, y = Height - 12 * q0;
            g.DrawLine(p, x0, y, x1, y); g.DrawLine(p, x0, y - 4 * q0, x0, y + 1); g.DrawLine(p, x1, y - 4 * q0, x1, y + 1);
            var sz = g.MeasureString("50 мм", f);
            g.DrawString("50 мм", f, Brushes.DimGray, x0 + (L - sz.Width) / 2, y - sz.Height - 2 * q0);
        }
    }
}

class MainForm : Form
{
    readonly RadioButton rbScale = new RadioButton { Text = "Масштаб конуса", Checked = true, AutoSize = true };
    readonly RadioButton rbFixed = new RadioButton { Text = "Постоянный размер", AutoSize = true };
    readonly TrackBar tbD = new TrackBar { Minimum = 50, Maximum = 314, Value = 100, TickFrequency = 25, SmallChange = 1, LargeChange = 10 };
    readonly TrackBar tbB = new TrackBar { Minimum = 0, Maximum = 12, Value = 4, TickFrequency = 1 };
    readonly Label lbD = new Label { AutoSize = true }, lbB = new Label { AutoSize = true }, lbF = new Label { AutoSize = true, ForeColor = Color.DimGray };
    readonly CheckBox cbRef = new CheckBox { Text = "Показать одну частоту 10 ГГц", Checked = true, AutoSize = true };
    readonly ListView lv = new ListView { View = View.Details, HeaderStyle = ColumnHeaderStyle.None, FullRowSelect = true, GridLines = false };
    readonly GeoView geo = new GeoView();
    readonly Plot pAz = new Plot { Title = "Сечение по азимуту", XLabel = "Азимут, град.", X0 = -180, X1 = 180, XT = new double[] { -180, -90, 0, 90, 180 } };
    readonly Plot pEl = new Plot { Title = "Сечение по углу места", XLabel = "Угол места, град.", X0 = -90, X1 = 90, XT = new double[] { -90, -45, 0, 45, 90 } };
    readonly Plot pRg = new Plot { Title = "Профиль по дальности", XLabel = "Дальность относительно цели, мм", X0 = -400, X1 = 400, XT = new double[] { -400, -200, 0, 200, 400 } };
    readonly double[] AZ, EL, DR; readonly double[][] dAz, dEl;
    readonly Timer debounce = new Timer { Interval = 40 };
    static readonly Color CUR = Color.FromArgb(15, 107, 92), REF = Color.FromArgb(154, 166, 180), WARN = Color.FromArgb(180, 83, 9);

    public MainForm()
    {
        Text = "Конус MIMO: шаг и полоса"; AutoScaleMode = AutoScaleMode.Dpi;
        float k = ScreenDpi() / 96f; Width = (int)(1280 * k); Height = (int)(800 * k); MinimumSize = new Size((int)(900 * k), (int)(600 * k));
        StartPosition = FormStartPosition.CenterScreen;
        try { Icon = Icon.ExtractAssociatedIcon(Application.ExecutablePath); } catch { }   // значок радара из exe
        Font = new Font("Segoe UI", 9.5f); BackColor = Color.FromArgb(243, 245, 247);
        AZ = Enumerable.Range(0, 721).Select(i => -180 + 0.5 * i).ToArray();
        EL = Enumerable.Range(0, 361).Select(i => -90 + 0.5 * i).ToArray();
        DR = Enumerable.Range(0, 401).Select(i => -400 + 2.0 * i).ToArray();
        dAz = AZ.Select(a => new[] { Math.Cos(a * Math.PI / 180), Math.Sin(a * Math.PI / 180), 0 }).ToArray();
        dEl = EL.Select(a => new[] { Math.Cos(a * Math.PI / 180), 0, Math.Sin(a * Math.PI / 180) }).ToArray();
        pAz.X = AZ; pEl.X = EL; pRg.X = DR;

        var left = new FlowLayoutPanel { Dock = DockStyle.Left, Width = (int)(330 * k), FlowDirection = FlowDirection.TopDown, WrapContents = false,
                                         Padding = new Padding(12), AutoScroll = true, BackColor = Color.White };
        Func<string, Label> H = s => new Label { Text = s.ToUpper(), AutoSize = true, ForeColor = Color.DimGray, Font = new Font("Segoe UI", 8.5f, FontStyle.Bold), Margin = new Padding(0, 10, 0, 4) };
        left.Controls.Add(H("Как менять шаг")); left.Controls.Add(rbScale); left.Controls.Add(rbFixed);
        left.Controls.Add(H("Шаг приёмников по дуге")); left.Controls.Add(lbD); tbD.Width = (int)(295 * k); left.Controls.Add(tbD);
        left.Controls.Add(H("Полоса сигнала")); left.Controls.Add(lbB); tbB.Width = (int)(295 * k); left.Controls.Add(tbB); left.Controls.Add(lbF);
        cbRef.Margin = new Padding(0, 10, 0, 0); left.Controls.Add(cbRef);
        left.Controls.Add(H("Что получилось"));
        lv.Width = (int)(295 * k); lv.Height = (int)(250 * k); lv.Columns.Add("", (int)(205 * k)); lv.Columns.Add("", (int)(80 * k), HorizontalAlignment.Right); left.Controls.Add(lv);

        var grid = new TableLayoutPanel { Dock = DockStyle.Fill, ColumnCount = 2, RowCount = 2, Padding = new Padding(8) };
        grid.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50)); grid.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50));
        grid.RowStyles.Add(new RowStyle(SizeType.Percent, 50)); grid.RowStyles.Add(new RowStyle(SizeType.Percent, 50));
        foreach (Control c in new Control[] { geo, pAz, pEl, pRg }) { c.Dock = DockStyle.Fill; c.Margin = new Padding(6); grid.Controls.Add(c); }
        Controls.Add(grid); Controls.Add(left);

        EventHandler ev = (s, e) => { debounce.Stop(); debounce.Start(); UpdateLabels(); };
        rbScale.CheckedChanged += ev; tbD.ValueChanged += ev; tbB.ValueChanged += ev; cbRef.CheckedChanged += ev;
        debounce.Tick += (s, e) => { debounce.Stop(); Compute(); };
        UpdateLabels(); Compute();
    }

    static float ScreenDpi() { using (var g = Graphics.FromHwnd(IntPtr.Zero)) return g.DpiX; }

    static string F(double v, int nd)
    {
        if (double.IsNaN(v) || double.IsInfinity(v)) return "–";
        return v.ToString("F" + nd, CultureInfo.InvariantCulture).Replace(".", ",").Replace("-", "−");
    }

    void UpdateLabels()
    {
        lbD.Text = "d = " + F(tbD.Value / 100.0, 2) + "λ";
        int b = tbB.Value * 5; var fs = Calc.Freqs(b);
        lbB.Text = b + " % от 10 ГГц";
        lbF.Text = b == 0 ? "одна частота 10 ГГц" : fs.Length + " частот, " + F(fs[0], 2) + "…" + F(fs[fs.Length - 1], 2) + " ГГц через 0,25 ГГц";
    }

    void Compute()
    {
        Cursor = Cursors.WaitCursor;
        double d = tbD.Value / 100.0; int b = tbB.Value * 5; bool scale = rbScale.Checked;
        var g = Calc.Geometry(d, scale); var fs = Calc.Freqs(b); var one = new[] { 10.0 };
        var az = Calc.Cut(g, fs, dAz); var el = Calc.Cut(g, fs, dEl); var rg = Calc.Range(fs, DR);
        bool showRef = cbRef.Checked && b > 0;
        double bwA, sA, bwE, sE, bwR, sR;
        Calc.Lobe(AZ, az, true, out bwA, out sA); Calc.Lobe(EL, el, false, out bwE, out sE); Calc.Lobe(DR, rg, false, out bwR, out sR);
        double cb = b > 0 ? Calc.C / (2 * 10.0 * b / 100) : double.NaN;

        geo.G = g; geo.Invalidate();
        pAz.Series.Clear(); pEl.Series.Clear(); pRg.Series.Clear(); pRg.VLines.Clear();
        if (showRef) { pAz.Series.Add(Tuple.Create(Calc.Cut(g, one, dAz), REF, true)); pEl.Series.Add(Tuple.Create(Calc.Cut(g, one, dEl), REF, true)); }
        pAz.Series.Add(Tuple.Create(az, CUR, false)); pEl.Series.Add(Tuple.Create(el, CUR, false)); pRg.Series.Add(Tuple.Create(rg, CUR, false));
        if (b > 0) { pRg.VLines.Add(Tuple.Create(-cb / 2, WARN)); pRg.VLines.Add(Tuple.Create(cb / 2, WARN)); }
        foreach (var pl in new[] { pAz, pEl })
        {
            pl.Legend.Clear(); pl.Legend.Add(Tuple.Create("выбранная полоса", CUR, false));
            if (showRef) pl.Legend.Add(Tuple.Create("одна частота", REF, true));
        }
        pRg.Legend.Clear(); pRg.Legend.Add(Tuple.Create("выбранная полоса", CUR, false));
        if (b > 0) pRg.Legend.Add(Tuple.Create("ширина c/(2B)", WARN, true));
        pAz.Invalidate(); pEl.Invalidate(); pRg.Invalidate();

        lv.BeginUpdate(); lv.Items.Clear();
        Action<string, string, bool> row = (k, v, bad) => { var it = new ListViewItem(new[] { k, v }); if (bad) it.ForeColor = WARN; lv.Items.Add(it); };
        row("Приёмники + передатчики", g.Rx.Count + " + " + g.Tx.Count, false);
        row("Виртуальных центров", (g.Rx.Count * g.Tx.Count).ToString(), false);
        row("Ø нижнего кольца, мм", F(2 * g.R1, 0), false);
        row("Шаг между кольцами, λ", F(g.RingStep / Calc.LAM, 2), false);
        row("Луч по азимуту, град.", F(bwA, 1), false);
        row("Луч по углу места, град.", F(bwE, 1), false);
        row("УБЛ по азимуту, дБ", F(sA, 1), sA > -13);
        row("УБЛ по углу места, дБ", F(sE, 1), sE > -13);
        row("Разрешение по дальности, мм", b > 0 ? F(bwR, 0) : "нет", false);
        row("Оценка c/(2B), мм", b > 0 ? F(cb, 0) : "–", false);
        lv.EndUpdate();
        Cursor = Cursors.Default;
    }

    [System.Runtime.InteropServices.DllImport("user32.dll")] static extern bool SetProcessDPIAware();

    [STAThread]
    static void Main(string[] args)
    {
        try { SetProcessDPIAware(); } catch { }         // чёткий текст при масштабе экрана 125-200 %
        Application.EnableVisualStyles(); Application.SetCompatibleTextRenderingDefault(false);
        var form = new MainForm();
        if (args.Length >= 2 && args[0] == "--shot")    // проверка: сохранить снимок окна и выйти
        {
            form.Shown += (s, e) =>
            {
                if (args.Length >= 5)                        // --shot файл scale|fixed d B%
                {
                    form.rbFixed.Checked = args[2] == "fixed";
                    form.tbD.Value = (int)Math.Round(double.Parse(args[3], CultureInfo.InvariantCulture) * 100);
                    form.tbB.Value = int.Parse(args[4]) / 5;
                    form.UpdateLabels(); form.Compute();
                }
                Application.DoEvents();
                using (var bmp = new Bitmap(form.Width, form.Height)) { form.DrawToBitmap(bmp, new Rectangle(0, 0, form.Width, form.Height)); bmp.Save(args[1]); }
                form.Close();
            };
        }
        Application.Run(form);
    }
}
