use probebwa::align::sensitive::{align_fixed_band, align_unbanded, SensitiveParameters};

#[test]
fn native_unbanded_sensitive_oracle() {
    check_sensitive(
        &std::env::var("STAMPY_SENSITIVE_FIXTURE")
            .map(|path| std::fs::read_to_string(path).unwrap())
            .unwrap_or_else(|_| include_str!("fixtures/stampy_sensitive.tsv").to_owned()),
    );
}

#[test]
fn native_fixed_band_sensitive_oracle() {
    check_sensitive(include_str!("fixtures/stampy_sensitive_band.tsv"));
}

fn check_sensitive(fixtures: &str) {
    let mut cases = 0;
    let mut differences = Vec::new();
    for line in fixtures
        .lines()
        .filter(|s| !s.starts_with('#') && !s.is_empty())
    {
        let f: Vec<_> = line.split('\t').collect();
        let quality: Vec<u8> = f[2].split(',').map(|s| s.parse().unwrap()).collect();
        let params = SensitiveParameters {
            insertion_prior: f[4].parse().unwrap(),
            gap_open: f[5].parse().unwrap(),
            gap_extend: f[6].parse().unwrap(),
            start_mean: f[7].parse().unwrap(),
            start_sd: f[8].parse().unwrap(),
        };
        let band: i32 = f.get(15).unwrap_or(&"-1").parse().unwrap();
        let got = if band < 0 {
            align_unbanded(f[1].as_bytes(), &quality, f[3].as_bytes(), params)
        } else {
            align_fixed_band(f[1].as_bytes(), &quality, f[3].as_bytes(), params, band)
        }
        .unwrap();
        let actual = [
            got.score.to_string(),
            String::from_utf8(got.primary.read).unwrap(),
            String::from_utf8(got.primary.reference).unwrap(),
            got.primary.offset.to_string(),
            String::from_utf8(got.alternate.read).unwrap(),
            String::from_utf8(got.alternate.reference).unwrap(),
        ];
        if actual.iter().map(String::as_str).collect::<Vec<_>>() != f[9..15] {
            differences.push(format!(
                "case {} expected {:?}, got {:?}",
                f[0],
                &f[9..15],
                actual
            ));
        }
        cases += 1;
    }
    assert!(cases >= 305);
    assert!(
        differences.is_empty(),
        "{}/{} cases differ:\n{}",
        differences.len(),
        cases,
        differences
            .iter()
            .take(12)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn native_coarse_oracle() {
    for line in include_str!("fixtures/stampy_coarse.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let qualities: Vec<u8> = f[2].split(',').map(|x| x.parse().unwrap()).collect();
        let got = probebwa::align::stampy_band::coarse_align(
            f[1].as_bytes(),
            &qualities,
            f[3].as_bytes(),
            f[4].parse().unwrap(),
            f[5].parse().unwrap(),
            f[6].parse().unwrap(),
        )
        .unwrap();
        assert_eq!(got.score.to_string(), f[7], "score case {}", f[0]);
        assert_eq!(got.read, f[8].as_bytes(), "read path case {}", f[0]);
        assert_eq!(
            got.reference,
            f[9].as_bytes(),
            "reference path case {}",
            f[0]
        );
        assert_eq!(got.offset.to_string(), f[10], "offset case {}", f[0]);
    }
}

#[test]
fn native_window_oracle() {
    let fixtures = std::env::var("STAMPY_WINDOW_FIXTURE")
        .map(|path| std::fs::read_to_string(path).unwrap())
        .unwrap_or_else(|_| decompress(include_bytes!("fixtures/stampy_windows.tsv.gz")));
    assert!(fixtures.lines().filter(|s| !s.starts_with('#')).count() >= 2630);
    for line in fixtures.lines().filter(|s| !s.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        let qualities: Vec<u8> = f[2].split(',').map(|x| x.parse().unwrap()).collect();
        let got = probebwa::align::stampy_band::find_window(
            f[1].as_bytes(),
            &qualities,
            f[3].as_bytes(),
            f[4].parse().unwrap(),
            f[5].parse().unwrap(),
            f[6].parse().unwrap(),
            f[7].parse().unwrap(),
        )
        .unwrap_or_else(|| panic!("no window event {}", f[0]));
        assert_eq!(
            (got.start, got.length, got.width),
            (
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
                f[10].parse().unwrap()
            ),
            "window event {}",
            f[0]
        );
    }
}

#[test]
fn native_sensitive_window_and_paths_oracle() {
    let fixtures = std::env::var("STAMPY_FULL_SENSITIVE_FIXTURE")
        .map(|path| std::fs::read_to_string(path).unwrap())
        .unwrap_or_else(|_| decompress(include_bytes!("fixtures/stampy_sensitive_full.tsv.gz")));
    let mut differences = Vec::new();
    let mut cases = 0;
    assert!(fixtures.lines().filter(|s| !s.starts_with('#')).count() >= 2630);
    for line in fixtures.lines().filter(|s| !s.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        let qualities: Vec<u8> = f[2].split(',').map(|x| x.parse().unwrap()).collect();
        let params = SensitiveParameters {
            insertion_prior: f[4].parse().unwrap(),
            gap_open: f[5].parse().unwrap(),
            gap_extend: f[6].parse().unwrap(),
            start_mean: f[7].parse().unwrap(),
            start_sd: f[8].parse().unwrap(),
        };
        let result = probebwa::align::sensitive::align_sensitive(
            f[1].as_bytes(),
            &qualities,
            f[3].as_bytes(),
            params,
            f[10].parse().unwrap(),
            f[9].parse().unwrap(),
        );
        cases += 1;
        let Some(got) = result else {
            differences.push(format!("no result event {}", f[0]));
            continue;
        };
        let offset = got.primary.offset;
        let (primary, alternate) = if f[11] == "3" || f[11] == "5" {
            (got.alternate, got.primary)
        } else {
            (got.primary, got.alternate)
        };
        let actual = [
            got.score.to_string(),
            String::from_utf8(primary.read).unwrap(),
            String::from_utf8(primary.reference).unwrap(),
            offset.to_string(),
            String::from_utf8(alternate.read).unwrap(),
            String::from_utf8(alternate.reference).unwrap(),
        ];
        if actual[..f.len() - 12]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            != f[12..]
        {
            differences.push(format!(
                "event {} expected {:?}, got {:?}",
                f[0],
                &f[12..],
                actual
            ));
        }
    }
    assert!(cases >= 263);
    assert!(
        differences.is_empty(),
        "{}/{} cases differ:\n{}",
        differences.len(),
        cases,
        differences
            .iter()
            .take(8)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

fn decompress(bytes: &[u8]) -> String {
    use std::io::Read;
    let mut result = String::new();
    flate2::read::GzDecoder::new(bytes)
        .read_to_string(&mut result)
        .unwrap();
    result
}

#[test]
fn native_baq_and_alignment_confidence() {
    let fixtures = decompress(include_bytes!("fixtures/stampy_confidence.tsv.gz"));
    let mut errors = Vec::new();
    let mut cases = 0;
    for line in fixtures.lines().filter(|s| !s.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        let read = f[1].as_bytes();
        let reference = f[3].as_bytes();
        let quals: Vec<u8> = f[2].split(',').map(|s| s.parse().unwrap()).collect();
        let mut p = SensitiveParameters {
            insertion_prior: f[4].parse().unwrap(),
            gap_open: f[5].parse().unwrap(),
            gap_extend: f[6].parse().unwrap(),
            start_mean: f[7].parse().unwrap(),
            start_sd: f[8].parse().unwrap(),
        };
        let window = probebwa::align::stampy_band::find_window(
            read,
            &quals,
            reference,
            p.gap_open as i32,
            p.gap_extend as i32,
            p.insertion_prior as i32,
            f[10]
                .parse::<i32>()
                .unwrap()
                .max(reference.len() as i32 - read.len() as i32),
        )
        .unwrap();
        p.start_mean -= window.start as i32;
        let reference = &reference[window.start..window.start + window.length];
        let a = align_fixed_band(read, &quals, reference, p, window.width).unwrap();
        let [first, second] = probebwa::align::confidence::confidence(
            read,
            &quals,
            reference,
            p,
            Some(window.width),
            [&a.primary, &a.alternate],
        )
        .unwrap();
        let probability = first.probability;
        let (primary, alternate) = if f[11] == "5" {
            (second, first)
        } else {
            (first, second)
        };
        let strings = [
            primary.quality_string,
            alternate.quality_string,
            primary.baq_offsets,
            alternate.baq_offsets,
        ];
        for (i, got) in strings.iter().enumerate() {
            if got != f[18 + i].as_bytes() {
                errors.push(format!(
                    "event {} field {} expected {} got {}",
                    f[0],
                    i,
                    f[18 + i],
                    String::from_utf8_lossy(got)
                ));
            }
        }
        let expected: f64 = f[22].parse().unwrap();
        if (probability - expected).abs() > 1e-10 {
            errors.push(format!(
                "event {} probability {} expected {}",
                f[0], probability, expected
            ));
        }
        cases += 1;
    }
    assert_eq!(cases, 209);
    assert!(
        errors.is_empty(),
        "{} discrepancies:\n{}",
        errors.len(),
        errors
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}
