## 0.10.0 / 2016-01-05

  * URI encode sitemap URLs (#85)
  * Do not include 'posts' collection twice (#92)
  * Fix GitHub Pages tests to test just the Jekyll version (#87)
  * Allow HTML files to end with `.xhtml` or `.htm` (#93)
  * Simplify whitespace regex for stripping whitespace (#96)

## 0.9.0 / 2015-09-21

  * Test against Jekyll 2, 3, and the GitHub Pages version. (#83)

## 0.8.1 / 2015-03-11

  * Do not assume all pages have changed (#35)
  * Remove duplicated range from regex (#73)

## 0.8.0 / 2015-02-03

  * Call each page `page` in pages loop in `sitemap.xml` for clarity (#64)
  * Remove `changefreq` (#34)
  * Remove `priority` (#33)
  * Don't strip 'index.html' when there is more to filename (#68)

## 0.7.0 / 2014-12-07

## Minor Enhancements

  * Make `site.baseurl` support more robust (#59)
  * Add `site.baseurl` to base site URL construction (#50)

## Development Fixes

  * Remove unnecessary spaces and escaping in README (#58)

## 0.6.3 / 2014-11-11

## Bug Fixes

  * Be backwards-compatible when `Site#in_source_dir` and `Site#in_dest_dir` don't exist (#57)

## 0.6.2 / 2014-11-08

### Bug Fixes

  * Don't attempt to read the sitemap upon page creation. (#52)
  * Use new secure methods to build source & dest paths. (#53)

## 0.6.1 / 2014-10-17

### Minor Enhancements

  * Strip excess whitespace (#40)

### Bug Fixes

  * Add UTC offset to `<lastmod>` to handle non-UTC timezones (#49)

### Development Fixes

  * Adding information about exclusion flag (#45)

## 0.6.0 / 2014-09-05

### Minor Enhancements

  * Include custom collections in the sitemap. (#30)
  * Use `post.last_modified_at` for post `<lastmod>` if available (#37)

## 0.5.1 / 2014-07-31

### Bug Fixes

  * Explicitly set sitemap layout to `nil` to avoid warning (#32)

## 0.5.0 / 2014-06-02

### Minor Enhancements

  * Allow users to exclude a page/post from the sitemap (#11)

## 0.4.1 / 2014-05-10

### Bug Fixes

  * Force sitemap layout to be `nil` (#16)
  * Correct seconds in timestamp for static files (#24)

# Development Fixes

  * Upgrade to Rspec 3.0 and use `be_truthy` (#24)

## 0.4.0 / 2014-05-06

### Major Enhancements

  * Support Jekyll 2.0 (#12)

## 0.3.0 / 2014-05-05

### Minor Enhancements

  * Generate sitemap using html_pages (#10)

### Bug Fixes

  * Remove stray sitemap.xsl from template (#8)

### Development Fixes

  * Added travis (#6)
  * Better timezone support (#7)

## 0.2.0 / 2014-03-24

  * Loosen Jekyll requirement (#4)

## 0.1.0 / 2014-03-15

  * Birthday!
