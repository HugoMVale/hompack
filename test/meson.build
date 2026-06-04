tests = [
	'TEST_F.f',
	'TEST_P.f',
	'TEST_S.f',
]

foreach t : tests
  test(
    t,
    executable(
      'test-@0@'.format(t.replace('.f', '')),
      sources: [t],
      link_with: hompack_lib,
    ),
    suite: 'unit',
    workdir: meson.project_source_root(),
  )
endforeach
