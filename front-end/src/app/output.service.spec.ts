import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { DomSanitizer } from '@angular/platform-browser';
import { TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';

import { BackendSocketClient } from './backend-socket.client';
import { Block } from './block.interface';
import { MockBackendSocketClient } from './mock-backend-socket.client';
import { MockUserIdService } from './mock-user-id.service';
import { OutputService } from './output.service';
import { UserIdService } from './user-id.service';

describe('OutputService', () => {
  let service: OutputService;
  let snackBar: MatSnackBar;
  let sanitizer: DomSanitizer;
  const clientMock = jasmine.createSpyObj('backend', ['listen', 'sendRequest']);

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        HttpClientTestingModule,
        MatSnackBarModule,
      ], 
      providers: [
        { provide: UserIdService, useClass: MockUserIdService },
        { provide: BackendSocketClient, useValue: clientMock}
      ],
    });
    snackBar = TestBed.inject(MatSnackBar);
    sanitizer = TestBed.inject(DomSanitizer);
  });

  it('should be created', () => {
    service = TestBed.inject(OutputService);
    expect(service).toBeTruthy();
  });

  describe('outputs', () => {
    beforeEach(()=>{
      service = TestBed.inject(OutputService);
    });

    it('should initially be empty', async () => {
      const outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('resetOutputs', () => {
    it('should set the outputs list to be empty', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(0);
      cb('{"text": "Dummy text"}');
      outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(1);
      expect(outputs[0].text).toBe('Dummy text');
      service.resetOutputs();
      outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('executeBlocks', () => {
    beforeEach(()=>{
      service = TestBed.inject(OutputService);
    });

    it('should open snack bar when userId is null', () => {
      const spy = spyOn(snackBar, 'open');
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      expect(spy).toHaveBeenCalledOnceWith('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    });

    it('should call backendSocketClient.sendRequest() with the appropriate message', () => {
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      backendSocketClient.sendRequest.calls.reset();
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {key: 'foo', text: 'bar', value: 42},
        ],
      };
      const blocks: Block[] = [];
      blocks.push(block);
      service.executeBlocks(blocks);
      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      expect(backendSocketClient.sendRequest).toHaveBeenCalledOnceWith(message);
    });
  });

  describe('Function passed to backendSocketClient.listen', () => {
    it('should add text to outputs array when it receives a valid text response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      cb('{"text": "Dummy text"}');
      outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(1);
      expect(outputs[0].text).toBe('Dummy text');
    });

    it('should add sanitised SafeUrl image and alt text to outputs array when it receives a valid image response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const spy = spyOn(sanitizer, 'bypassSecurityTrustUrl');
      cb('{"img": "Image text", "alttext": "Alt text"}');
      outputs = await firstValueFrom(service.outputs);
      expect(outputs.length).toBe(1);
      expect(spy).toHaveBeenCalledTimes(1);
      const expectedValue = sanitizer.bypassSecurityTrustUrl('data:image/png;base64,' + 'Image text');
      expect(outputs[0].img).toBe(expectedValue);
      expect(outputs[0].alttext).toBe('Alt text');
    });

    it('should not change the outputs array when it receives an end_connection response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      cb('{"end_connection": "end_connection"}');
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });

    it('should open snack bar with given message when response is a known error', async() => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const spy = spyOn(snackBar, 'open');
      cb('{"error": "error"}');
      expect(spy).toHaveBeenCalledOnceWith('error', 'Close', { duration: 5000 });
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });

    it('should open snack bar with generic message for invalid responses', async() => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
      const spy = spyOn(snackBar, 'open');
      cb('{}');
      expect(spy).toHaveBeenCalledOnceWith('Received a bad response, please refresh the page and try again', 'Close', { duration: 5000 });
      outputs = await firstValueFrom(service.outputs);
      expect(outputs).toEqual([]);
    });
  });

  describe('executingBlocks', () => {
    it('should be false by default', async () => {
      service = TestBed.inject(OutputService);
      const executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
    });

    it('should be set to true by executeBlocks function', async () => {
      service = TestBed.inject(OutputService);
      let executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeTrue();
    });

    it('should be set to false after an end_connection response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeTrue();
      cb('{"end_connection": "end_connection"}');
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
    });

    it('should be set to false after a known error response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeTrue();
      cb('{"error": "error"}');
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
    });

    it('should be set to false after an invalid response', async () => {
      const backendSocketClient: jasmine.SpyObj<BackendSocketClient> = TestBed.inject(BackendSocketClient) as jasmine.SpyObj<BackendSocketClient>;
      let cb: (res: string) => void = () => { throw Error('Should have been changed'); };
      backendSocketClient.listen.and.callFake(func => {
        cb = func;
      });
      service = TestBed.inject(OutputService);
      expect(backendSocketClient.listen).toHaveBeenCalled();
      let executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeTrue();
      cb('{}');
      executingBlocks = await firstValueFrom(service.executingBlocks);
      expect(executingBlocks).toBeFalse();
    });
  });
});
