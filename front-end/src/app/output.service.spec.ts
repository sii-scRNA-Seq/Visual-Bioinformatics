import { TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { OutputService } from './output.service';

describe('OutputService', () => {
  let service: OutputService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule]
    });
    service = TestBed.inject(OutputService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
