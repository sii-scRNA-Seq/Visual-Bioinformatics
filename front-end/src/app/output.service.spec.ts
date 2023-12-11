import { TestBed, fakeAsync, tick } from '@angular/core/testing';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { OutputService } from './output.service';
import { first, firstValueFrom, from } from 'rxjs';
import { Block } from './block.service';

describe('OutputService', () => {
  let service: OutputService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [HttpClientTestingModule],
    });
    service = TestBed.inject(OutputService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('Outputs', () => {
    it('should initially be empty', async () => {

      const outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
      
    });
  });

  describe('Execute Block', () => {
    it('should add a response to outputs array', fakeAsync(() => {
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      const mockHttp = TestBed.inject(HttpTestingController);
      service.executeBlock(block).then(async () => {
        const outputs = await firstValueFrom(service.outputs);
        expect(outputs).toEqual([{text: 'Hello', other: 'World'}]);
      });
      tick();
      const req = mockHttp.expectOne('http://127.0.0.1:5000/loaddata/');
      expect(req.request.method).toBe('GET');
      req.flush({text: 'Hello', other: 'World'});
      mockHttp.verify();  

    }));
  }); 
});
